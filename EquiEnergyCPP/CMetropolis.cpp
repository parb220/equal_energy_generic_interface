#include "CMetropolis.h"
#include <cmath>
#include "dw_rand.h"
#include <fstream>

using namespace std; 

void CMetropolis::BlockAdaptive(const CSampleIDWeight &adaptive_start_point, const vector<TDenseMatrix> &B, double mid, size_t period, size_t max_period)
// Makes block-wise random-walk Metropolis draws, adjusting the block scale until the target acceptance rate 
// is hit.
//
// 	adaptive_start_point:	n*1 column vector, initial value
// 	B:			k matrices, each of which is n*b(i) consisting of the direction vectors of the i-th block
// 	mid:			target acceptance ratio
// 	period:			initial number of draws to make before recomputing scale. Increases by a factor of 2
// 				when target ratio has been reached
// 	max_period:		final number of draws to make before recomputing scale. Should be a power of 2 times
// 				period
//
// For each draw, the block i jumping kernel is as follows
// 	(1) draw z from the b(i)-dimensional standard normal distribution
// 	(2) let x = scale(i)*B{i}*z
// 	(3) if y is the current point, then y+x is the proposal point
//
// Usually b(1)+...+b(k) = n and the matrix[B{1} ... B{k}] has orthogonal columsn, though this is not checked	
{
	size_t k = B.size(); 	// number of blocks

	unsigned int n_draws =0;
	vector<unsigned int> n_accepted(k,0); 
	vector<unsigned int> begin_draws(k,0); 
	vector<unsigned int> begin_jumps(k,0); 
	vector<unsigned int> periods(k,period); 
	vector<unsigned int> end_draws(k,period); 
	
	TDenseVector previous_ratio(k,0.0); 
	TDenseVector scale(k,1.0); 
	TDenseVector best_scale(k,1.0); 
	TDenseVector low_scale(k,-1.0); 
	TDenseVector low_jump_ratio(k,-1.0); 
	TDenseVector high_scale(k,-1.0); 
	TDenseVector high_jump_ratio(k,-1.0); 

	double log_mid = log(mid); 
	double lower_bound = exp(log_mid/0.2); 	
	double upper_bound = exp(log_mid/5.0);

	vector<size_t> b(k,0); // b(i): size of the i-th block
	for (unsigned int i=0; i<k; i++)
		b[i] = B[i].cols;  

	// Beginning adaptive burn-in 
	CSampleIDWeight x, y; 
	y=adaptive_start_point; 
	double log_previous = model->log_posterior_function(y), log_current, new_scale, diff;
	bool done = false; 
	unsigned int check = period; 
	while (!done)
	{
		// draw metropolis blocks
		for (unsigned int i=0; i<k; i++)
		{
			x.data = y.data + scale[i]*(B[i]*RandomNormalVector(b[i]));
			log_current = model->log_posterior_function(x); 

			if (log_current-log_previous >= log(dw_uniform_rnd()) )
			{
				y= x; 
				log_previous = log_current; 
				n_accepted[i] ++; 
			}
		}
		n_draws ++; 

		// rescale
		if (n_draws == check)
		{
			done = true; 
			for (unsigned int i=0; i<k; i++)
			{
				// jump ratio
				previous_ratio[i] = (double)(n_accepted[i]-begin_jumps[i])/(double)(n_draws-begin_draws[i]); 
				// recompute scale?
				if (end_draws[i] <= n_draws)
				{
					// set new low or high bounds
					if (previous_ratio[i] < mid)
					{
						low_scale[i] = scale[i]; 
						low_jump_ratio[i] = previous_ratio[i]; 
					}
					else 
					{
						high_scale[i] = scale[i]; 
						high_jump_ratio[i] = previous_ratio[i]; 
					}

					// new scale and best scale
					if (low_jump_ratio[i] < 0.0)
					{
						best_scale[i] = scale[i]; 
						if (previous_ratio[i] > upper_bound)
							new_scale = 5.0*high_scale[i]; 
						else 
							new_scale = (log_mid/log(previous_ratio(i)))*high_scale[i];
					}
					else if (high_jump_ratio[i] < 0.0)
					{
						best_scale[i] = scale[i]; 
						if (previous_ratio[i] < lower_bound)
							new_scale = 0.2*low_scale[i]; 
						else 
							new_scale=(log_mid/log(previous_ratio(i)))*low_scale[i];
					}
					else 
					{
						diff = high_jump_ratio[i]-low_jump_ratio[i]; 
						if (diff > 1.0e-6)
							new_scale=((mid-low_jump_ratio[i])*low_scale[i] + (high_jump_ratio[i]-mid)*high_scale[i])/diff;
						else 
							new_scale=0.5*(low_scale[i] + high_scale[i]);
						best_scale[i] = new_scale; 
						periods[i] = 2* periods[i]; 
						low_jump_ratio[i] = -1.0; 
						high_jump_ratio[i] = -1.0; 
					}
			
					// reset adaptive counts and scale
					begin_jumps[i] = n_accepted[i]; 
					begin_draws[i] = n_draws; 
					end_draws[i] = n_draws+periods[i]; 
					scale[i] = new_scale; 
				}

				// not done?
				if (periods[i] <= max_period)
					done = false; 
			}
			check += period; 
		}
	}
	
	blocks.resize(k); 
	for (unsigned int i=0; i<k; i++)
		blocks[i] = B[i]*scale[i]; 
}

bool CMetropolis:: BlockRandomWalkMetropolis(double &log_posterior_y, CSampleIDWeight &y, const CSampleIDWeight &initial_v)
// Given
//	initial_v:	initial value, vector
// Results:
// 	y:	new vector according to random-walk metropolis
// 	log_posterior_y:	log_posterior(y)
// Returns:
// 	true:	if y is a new sample (different than x)
// 	false:	if y is the same as x
// Involving:
// 	blocks:	k matrices, each of which is n*b(i) consisting of the direction vectors of the i-th block
// 	model:	contains target distribution	
{
	size_t k=blocks.size(); 	// number of blocks	
	vector<size_t> b; 	// size of each block
	for (unsigned int i=0; i<k; i++)
		b[i] = blocks[i].cols; 

	CSampleIDWeight x = initial_v;  
	double log_previous = model->log_posterior_function(x), log_current; 
	bool if_new_sample = false; 
	for (unsigned int i=0; i<k; i++)
	{
		y.data = x.data+blocks[i]*RandomNormalVector(b[i]); 
		log_current = model->log_posterior_function(y); 
		if (log_current - log_previous >= log(dw_uniform_rnd()) )
		{
			x = y; 
			log_previous = log_current; 
			if (!if_new_sample)
				if_new_sample = true; 
		}
	}
	
	log_posterior_y = log_previous; 
	return if_new_sample; 
}

void CMetropolis::FourPassAdaptive(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin)
// first pass
// 	identity variance matrix / n blocks
// second pass
// 	diagonal variance matrix / 1 block
// third pass 
// 	simulate from second pass block structure
// 	variance matrix estimated from simulation / n blocks
// forth pass
// 	variance matrix estimated from simulation with scaled columns / 1 block
{
	CSampleIDWeight x=adaptive_start_point; 

	// first pass
	size_t n_blocks = x.data.dim; 
	vector<TDenseMatrix> B_matrix(n_blocks); 
	TDenseMatrix I_matrix = Identity(x.data.dim); 
	for (unsigned int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(I_matrix,i));
	BlockAdaptive(x, B_matrix, 0.25, period, max_period); 
	// blocks should now contain n matrixes, each of which is a n-by-1 matrix

	// second pass
	n_blocks = 1; 
	B_matrix.resize(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}(:,i) = blocks{i} <=> B[0][:,i] = blocks[i][:,0]
	for (unsigned int i=0; i<x.data.dim; i++)
		for (unsigned int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0);  		
	BlockAdaptive(x, B_matrix, 0.25, period, max_period); 

	// simulate
	// burn-in
	double log_posterior_y; 
	CSampleIDWeight y; 
	for (unsigned int t=0; t<burn_in; t++)
	{
		if (BlockRandomWalkMetropolis(log_posterior_y, y, x) )
			x = y; 
	}

	vector<CSampleIDWeight> Y_simulation(n_draws); 
	// simulate
	for (unsigned int i=0; i<n_draws; i++)
	{
		for (unsigned int t=0; t<thin; t++)
		{
			if (BlockRandomWalkMetropolis(log_posterior_y, y, x) )
                        	x = y; 
		}
		Y_simulation[i]=y; 
	} 

	// Compute variance and mean
	TDenseVector mean(y.data.dim,0.0);  
	TDenseMatrix variance(y.data.dim,y.data.dim,0.0), U_matrix, V_matrix, D_matrix; 
	TDenseVector d_vector; 
	for (unsigned int i=0; i<n_draws; i++)
	{
		mean = mean + Y_simulation[i].data; 
		variance = variance + Multiply(Y_simulation[i].data,Y_simulation[i].data);
	}
	
	mean = (1.0/(double)n_draws) * mean; 
	variance = (1.0/(double)n_draws)*(variance+Transpose(variance)); 
	variance = variance - Multiply(mean, mean); 
 
	SVD(U_matrix, d_vector, V_matrix, variance); 
	D_matrix = DiagonalMatrix(d_vector); 
	U_matrix = U_matrix *D_matrix; 

	// third-pass: n blocks
	n_blocks = x.data.dim; 
	B_matrix.resize(n_blocks); 
	for (unsigned int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(U_matrix,i)); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period);
	// blocks should contain n_blocks matrices, each of them being a x.dim-by-1 matrix
	
	// forth-pass: 1 block
	n_blocks = 1; 
	B_matrix.resize(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}{:,i} = blocks{i} <=> B[0][:,i] = blocks[i](:,0)
	for (unsigned int i=0; i<x.data.dim; i++)
		for (unsigned int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period); 
}

bool CMetropolis::ReadBlocks(const string &file_name)
{
	fstream iFile(file_name.c_str(), ios::in|ios::binary); 
	if (!iFile)
		return false; 
	size_t n_blocks; 
	iFile.read((char*)(&n_blocks),sizeof(size_t) ); 
	blocks.resize(n_blocks); 
	for (unsigned int i=0; i<n_blocks; i++)
	{
		size_t m, n; 
		iFile.read((char*)(&m),sizeof(size_t) ); 
		iFile.read((char*)(&n),sizeof(size_t) ); 
		blocks[i].Resize(m,n); 
		for (unsigned int j=0; j<n; j++)
			for (unsigned int k=0; k<m; k++)
			{
				double data; 
				iFile.read((char*)(&data),sizeof(double)); 
				blocks[i](k,j)=data; 
			}	
	}

	iFile.close(); 
	return true; 	
}

bool CMetropolis::WriteBlocks(const string &file_name)
{
	fstream oFile(file_name.c_str(), ios::out|ios::binary); 
	if (!oFile)
		return false;
	size_t n_blocks = blocks.size(); 
	oFile.write((char*)(&n_blocks), sizeof(size_t));	// number of blocks
	for (unsigned int i=0; i<n_blocks; i++)
	{
		size_t m = blocks[i].rows; 
		size_t n = blocks[i].cols; 
		oFile.write((char*)(&m), sizeof(size_t)); 
		oFile.write((char*)(&n), sizeof(size_t)); 
		for (unsigned int j=0; j<n; j++)
			for (unsigned int k=0; k<m; k++)
			{
				double data = blocks[i](k,j); 
				oFile.write((char*)(&data), sizeof(double)); 
			}
	}
	oFile.close(); 
	return true;
}

bool CMetropolis::AdaptiveBeforeSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &block_file_name)
// There are two passes before simulations
// first pass
// 	identity variance matrix / n blocks
// second pass
// 	diagonal variance matrix / 1 block
// Here file_name refers to a file that blocks are written into. 
// Simulation will first read blocks from the file and then simulate
{
	CSampleIDWeight x=adaptive_start_point; 

	// first pass
	size_t n_blocks = x.data.dim; 
	vector<TDenseMatrix> B_matrix(n_blocks); 
	TDenseMatrix I_matrix = Identity(x.data.dim); 
	for (unsigned int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(I_matrix,i));
	BlockAdaptive(x, B_matrix, 0.25, period, max_period); 
	// blocks should now contain n matrixes, each of which is a n-by-1 matrix

	// second pass
	n_blocks = 1; 
	B_matrix.resize(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}(:,i) = blocks{i} <=> B[0][:,i] = blocks[i][:,0]
	for (unsigned int i=0; i<x.data.dim; i++)
		for (unsigned int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0);  		
	BlockAdaptive(x, B_matrix, 0.25, period, max_period); 

	return WriteBlocks(block_file_name); 
}

bool CMetropolis::AdaptiveAfterSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &sample_file_name, const string &block_file_name) 
// There are two more passes after Simulation
// third pass 
// 	variance matrix estimated from simulation / n blocks
// forth pass
// 	variance matrix estimated from simulation with scaled columns / 1 block
// Here filename refers to the file where samples are stored to. 
// So we first load these samples, and then calculate mean and variance. 
// finally blocks will be written into a file so that simulation can start from there
{	
	vector<CSampleIDWeight> Y; 
	if (!LoadSampleFromFile(sample_file_name, Y) )
		return false; 		

	// Compute variance and mean
	TDenseVector mean(Y[0].data.dim,0.0);  
	TDenseMatrix variance(Y[0].data.dim,Y[0].data.dim,0.0), U_matrix, V_matrix, D_matrix; 
	TDenseVector d_vector; 
	for (unsigned int i=0; i<Y.size(); i++)
	{
		mean = mean + Y[i].data; 
		variance = variance + Multiply(Y[i].data,Y[i].data);
	}
	
	mean = (1.0/(double)Y.size()) * mean; 
	variance = (1.0/(double)Y.size())*(variance+Transpose(variance)); 
	variance = variance - Multiply(mean, mean); 
 
	SVD(U_matrix, d_vector, V_matrix, variance); 
	D_matrix = DiagonalMatrix(d_vector); 
	U_matrix = U_matrix *D_matrix; 

	// third-pass: n blocks
	CSampleIDWeight x=adaptive_start_point; 
	size_t n_blocks = x.data.dim; 
	vector<TDenseMatrix> B_matrix(n_blocks); 
	for (unsigned int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(U_matrix,i)); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period);
	
	// forth-pass: 1 block
	n_blocks = 1; 
	B_matrix.resize(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}{:,i} = blocks{i} <=> B[0][:,i] = blocks[i](:,0)
	for (unsigned int i=0; i<x.data.dim; i++)
		for (unsigned int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period); 
	return WriteBlocks(block_file_name); 
}

