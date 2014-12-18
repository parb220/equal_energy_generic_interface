#include <cmath>
#include <fstream>
#include "dw_rand.h"
#include "CSampleIDWeight.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"

using namespace std; 

void CMetropolis::BlockAdaptive(const CSampleIDWeight &adaptive_start_point, const vector<TDenseMatrix> &B, double mid, size_t period, size_t max_period, bool if_eejump)
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
	blocks = vector<TDenseMatrix>(k);  
	
	int n_draws =0;
	vector<int> n_accepted(k,0); 
	vector<int> begin_draws(k,0); 
	vector<int> begin_jumps(k,0); 
	vector<int> periods(k,period); 
	vector<int> end_draws(k,period); 
	
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
	for (int i=0; i<k; i++)
		b[i] = B[i].cols;  

	// Beginning adaptive burn-in 
	CSampleIDWeight y=adaptive_start_point, x; 
	double log_previous = model->log_posterior_function(y), log_current, new_scale, diff;
	bool done = false; 
	int check = period; 
	while (!done)
	{
		if (if_eejump && model->MakeEquiEnergyJump(x,y))
		{
			y = x; 
			log_previous = model->log_posterior_function(y); 
		}
		// draw metropolis blocks
		for (int i=0; i<k; i++)
		{
			x.data = y.data + scale[i]*(B[i]*RandomNormalVector(b[i]));
			x.DataChanged(); 
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
			for (int i=0; i<k; i++)
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
	
	for (int i=0; i<k; i++)
		blocks[i] = B[i]*best_scale[i]; 
}

bool CMetropolis:: BlockRandomWalkMetropolis(double &log_posterior_y, CSampleIDWeight &y, const CSampleIDWeight &initial_v, size_t thin)
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
	vector<size_t> b(k); 	// size of each block
	for (int i=0; i<k; i++)
		b[i] = blocks[i].cols; 

	CSampleIDWeight x = initial_v;  
	double log_previous = model->log_posterior_function(x), log_current; 
	bool if_new_sample = false;  
	for (int i_thin=0; i_thin<thin; i_thin++)
	{
		for (int i=0; i<k; i++)
		{
			y.data = x.data+blocks[i]*RandomNormalVector(b[i]); 
			y.DataChanged(); 
			log_current = model->log_posterior_function(y); 
			if (log_current - log_previous >= log(dw_uniform_rnd()) )
			{	
				x = y; 
				log_previous = log_current; 
				if (!if_new_sample)
					if_new_sample = true; 
			}
		}
	}
	
	log_posterior_y = log_previous; 
	y = x; 
	return if_new_sample; 
}

bool CMetropolis::FourPassAdaptive_StartWithoutSampleFile(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, const string &block_file_name, bool if_eejump)
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
	for (int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(I_matrix,i));
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump); 
	// blocks should now contain n matrixes, each of which is a n-by-1 matrix

	// second pass
	n_blocks = 1; 
	B_matrix = vector<TDenseMatrix>(n_blocks);
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}(:,i) = blocks{i} <=> B[0][:,i] = blocks[i][:,0]
	for (int i=0; i<x.data.dim; i++)
		for (int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0);  		
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump); 

	// simulate
	// burn-in
	double log_posterior_y; 
	CSampleIDWeight y; 
	for (int t=0; t<burn_in; t++)
	{
		if (BlockRandomWalkMetropolis(log_posterior_y, y, x) )
			x = y; 
	}

	vector<CSampleIDWeight> Y(n_draws); 
	// simulate
	for (int i=0; i<n_draws; i++)
	{
		for (int t=0; t<thin; t++)
		{
			if (BlockRandomWalkMetropolis(log_posterior_y, y, x) )
                        	x = y; 
		}
		Y[i]=x; 
	} 

	// Compute variance and mean
	TDenseVector mean(y.data.dim,0.0);  
	TDenseMatrix variance(y.data.dim,y.data.dim,0.0), U_matrix, V_matrix, D_matrix; 
	TDenseVector d_vector; 
	for (int i=0; i<n_draws; i++)
	{
		mean = mean + Y[i].data; 
		variance = variance + Multiply(Y[i].data,Y[i].data);
	}
	
	mean = (1.0/(double)n_draws) * mean; 
	variance = (1.0/(double)n_draws)*variance - Multiply(mean, mean);
	variance = 0.5*(variance+Transpose(variance));  
 
	SVD(U_matrix, d_vector, V_matrix, variance); 
	for (int i=0; i<d_vector.dim; i++)
		d_vector[i] = sqrt(d_vector[i]); 
	D_matrix = DiagonalMatrix(d_vector); 
	U_matrix = U_matrix *D_matrix; 

	// third-pass: n blocks
	n_blocks = x.data.dim; 
	B_matrix = vector<TDenseMatrix>(n_blocks); 
	for (int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(U_matrix,i)); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump);
	// blocks should contain n_blocks matrices, each of them being a x.dim-by-1 matrix
	
	// forth-pass: 1 block
	n_blocks = 1; 
	B_matrix = vector<TDenseMatrix>(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}{:,i} = blocks{i} <=> B[0][:,i] = blocks[i](:,0)
	for (int i=0; i<x.data.dim; i++)
		for (int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump); 
	if (block_file_name.empty())
		return false; 
	return WriteBlocks(block_file_name); 
}

bool CMetropolis::ReadBlocks(const string &file_name)
{
	fstream iFile(file_name.c_str(), ios::in|ios::binary); 
	if (!iFile)
		return false; 
	size_t n_blocks; 
	iFile.read((char*)(&n_blocks),sizeof(size_t) ); 
	blocks = vector<TDenseMatrix>(n_blocks); 
	for (int i=0; i<n_blocks; i++)
	{
		size_t m, n; 
		iFile.read((char*)(&m),sizeof(size_t) ); 
		iFile.read((char*)(&n),sizeof(size_t) ); 
		blocks[i]=TDenseMatrix(m,n); 
		for (int j=0; j<n; j++)
			for (int k=0; k<m; k++)
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
	for (int i=0; i<n_blocks; i++)
	{
		size_t m = blocks[i].rows; 
		size_t n = blocks[i].cols; 
		oFile.write((char*)(&m), sizeof(size_t)); 
		oFile.write((char*)(&n), sizeof(size_t)); 
		for (int j=0; j<n; j++)
			for (int k=0; k<m; k++)
			{
				double data = blocks[i](k,j); 
				oFile.write((char*)(&data), sizeof(double)); 
			}
	}
	oFile.close(); 
	return true;
}

bool CMetropolis::AggregateBlocksAndRemoveFiles(const vector<string> &read_file_name, const string &write_file_name)
{
	if (read_file_name.size() == 1)
		rename(read_file_name[0].c_str(), write_file_name.c_str());	
	else 
	{
		vector<TDenseMatrix> sum_blocks; 
		for (int i=0; i<(int)read_file_name.size(); i++)
		{
			if (!ReadBlocks(read_file_name[i]) )
				return false; 
			if (i == 0)
			{
				sum_blocks.resize(blocks.size()); 
				for (int j=0; j<(int)blocks.size(); j++)
					sum_blocks[j] = blocks[j]; 
			}
			else 
			{
				try { 
					for (int j=0; j<(int)(blocks.size()); j++)
						sum_blocks[j] = sum_blocks[j] + blocks[j]; 
				}
				catch(...) { return false; }
			}
			remove(read_file_name[i].c_str()); 
		}
		for (int j=0; j<(int)(sum_blocks.size()); j++)
		{
			sum_blocks[j] = sum_blocks[j] * (1.0/(double)read_file_name.size()); 
			blocks[j] = sum_blocks[j]; 
		}
		if (!WriteBlocks(write_file_name))
			return false; 
	}
	return true; 
}

bool CMetropolis::AdaptiveBeforeSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &block_file_name, bool if_eejump)
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
	for (int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(I_matrix,i));
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump); 
	// blocks should now contain n matrixes, each of which is a n-by-1 matrix

	// second pass
	n_blocks = 1; 
	B_matrix = vector<TDenseMatrix>(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}(:,i) = blocks{i} <=> B[0][:,i] = blocks[i][:,0]
	for (int i=0; i<x.data.dim; i++)
		for (int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0);  		
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump); 

	if (block_file_name.empty())
		return true; 
	else 
	{
		bool return_value = WriteBlocks(block_file_name); 
		return return_value; 
	}
}

bool CMetropolis::AdaptiveAfterSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &sample_file_name, const string &block_file_name, bool if_eejump) 
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
	for (int i=0; i<Y.size(); i++)
	{
		mean = mean + Y[i].data; 
		variance = variance + Multiply(Y[i].data,Y[i].data);
	}
	
	mean = (1.0/(double)Y.size()) * mean; 
	variance = (1.0/(double)Y.size())*variance - Multiply(mean, mean); 
	variance = 0.5 * (variance+Transpose(variance)); 
 
	SVD(U_matrix, d_vector, V_matrix, variance); 
	for (int i=0; i<d_vector.dim; i++)
		d_vector[i] = sqrt(d_vector[i]); 
	D_matrix = DiagonalMatrix(d_vector); 
	U_matrix = U_matrix *D_matrix; 

	// third-pass: n blocks
	CSampleIDWeight x=adaptive_start_point; 
	size_t n_blocks = x.data.dim; 
	vector<TDenseMatrix> B_matrix(n_blocks); 
	for (int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(U_matrix,i)); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump);
	
	// forth-pass: 1 block
	n_blocks = 1; 
	B_matrix = vector<TDenseMatrix>(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}{:,i} = blocks{i} <=> B[0][:,i] = blocks[i](:,0)
	for (int i=0; i<x.data.dim; i++)
		for (int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump); 

	if (block_file_name.empty())
		return true; 
	else 
	{
		bool return_value = WriteBlocks(block_file_name); 
		return return_value; 
	}
}

bool CMetropolis::FourPassAdaptive_StartWithSampleFile(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, const string &sample_file_name, const string &block_file_name, bool if_eejump)
// third pass (1)
// 	variance matrix estimated from the given file / n blocks
// fourth pass (1)
// 	variance matrix estimated from the given file / 1 block
// third pass (2)
// 	simulate from third pass (1) block structure
// 	variance matrix estimated from simulation / n blocks
// fourth pass (2)
// 	variance matrix estimated from simulation with scaled columns / 1 block
{
	vector<CSampleIDWeight> Y; 
	if (!LoadSampleFromFile(sample_file_name, Y) )
		return false; 		

	// Compute variance and mean
	TDenseVector mean(Y[0].data.dim,0.0);  
	TDenseMatrix variance(Y[0].data.dim,Y[0].data.dim,0.0), U_matrix, V_matrix, D_matrix; 
	TDenseVector d_vector; 
	for (int i=0; i<Y.size(); i++)
	{
		mean = mean + Y[i].data; 
		variance = variance + Multiply(Y[i].data,Y[i].data);
	}
	
	mean = (1.0/(double)Y.size()) * mean; 
	variance = (1.0/(double)Y.size())*variance - Multiply(mean, mean); 
	variance = 0.5*(variance+Transpose(variance)); 
 
	SVD(U_matrix, d_vector, V_matrix, variance); 
	for (int i=0; i<d_vector.dim; i++)
		d_vector[i] = sqrt(d_vector[i]); 
	D_matrix = DiagonalMatrix(d_vector); 
	U_matrix = U_matrix *D_matrix; 

	// third-pass (1): n blocks
	CSampleIDWeight x=adaptive_start_point; 
	size_t n_blocks = x.data.dim; 
	vector<TDenseMatrix> B_matrix(n_blocks); 
	for (int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(U_matrix,i)); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump);
	
	// forth-pass (1) : 1 block
	n_blocks = 1; 
	B_matrix = vector<TDenseMatrix>(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}{:,i} = blocks{i} <=> B[0][:,i] = blocks[i](:,0)
	for (int i=0; i<x.data.dim; i++)
		for (int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0); 

	// simulate
	// burn-in
	double log_posterior_y; 
	CSampleIDWeight y; 
	for (int t=0; t<burn_in; t++)
	{
		if (BlockRandomWalkMetropolis(log_posterior_y, y, x) )
			x = y; 
	}

	Y = vector<CSampleIDWeight>(n_draws); 
	// simulate
	for (int i=0; i<n_draws; i++)
	{
		for (int t=0; t<thin; t++)
		{
			if (BlockRandomWalkMetropolis(log_posterior_y, y, x) )
                        	x = y; 
		}
		Y[i]=x; 
	} 

	// Compute variance and mean
	mean.Zeros(); 
	variance.Zeros(); 
	for (int i=0; i<n_draws; i++)
	{
		mean = mean + Y[i].data; 
		variance = variance + Multiply(Y[i].data,Y[i].data);
	}
	
	mean = (1.0/(double)n_draws) * mean; 
	variance = (1.0/(double)n_draws)*variance - Multiply(mean, mean); 
	variance = 0.5*(variance+Transpose(variance)); 
 
	SVD(U_matrix, d_vector, V_matrix, variance); 
	for (int i=0; i<d_vector.dim; i++)
		d_vector[i] = sqrt(d_vector[i]); 
	D_matrix = DiagonalMatrix(d_vector); 
	U_matrix = U_matrix *D_matrix; 

	// third-pass (2): n blocks
	n_blocks = x.data.dim; 
	B_matrix = vector<TDenseMatrix>(n_blocks); 
	for (int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(U_matrix,i)); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump);
	// blocks should contain n_blocks matrices, each of them being a x.dim-by-1 matrix
	
	// forth-pass(2): 1 block
	n_blocks = 1; 
	B_matrix = vector<TDenseMatrix>(n_blocks); 
	B_matrix[0].Zeros(x.data.dim,x.data.dim); 
	// B{1}{:,i} = blocks{i} <=> B[0][:,i] = blocks[i](:,0)
	for (int i=0; i<x.data.dim; i++)
		for (int j=0; j<x.data.dim; j++)
			B_matrix[0](j,i) = blocks[i](j,0); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period, if_eejump); 
	if (block_file_name.empty())
		return false; 
	return WriteBlocks(block_file_name); 
}

bool CMetropolis::AdaptiveBeforeSimulation_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &block_file_name, bool if_eejump, const string &block_scheme_file_name)
// block_scheme contains the assignment of the dimensions into the blocks. 
//
// If there is only one block, then B_matrix only contains one identity matrix.
// If there are k blocks, the B_matrix contains k matrices. Each matrix is a n-by-bi matrix  where the each 
//                column corresponds to a dimension that has been assigned to the i-th block. Let e_j be a 
//                vector whose j-th element is 1 and all other elements are 0's. 
//                Then B_matrix[i] = [e_1(i), e_2(i), ...,], with 1(i) indicates the 1st dimension assigned to 
//                the i-th block. 
// Here block_file_name refers to the file that blocks are written into. 
// block_scheme_file_name refers to the file that contains the assignment of the dimensions into the blocks.
// Simulation will first read blocks from the file and then simulate
{
	CSampleIDWeight x=adaptive_start_point; 

	block_scheme = ReadBlockScheme(block_scheme_file_name); 
	if (block_scheme.empty())
		block_scheme.push_back(TIndex(0,x.data.dim-1)); 

	vector<TDenseMatrix> B_matrix(block_scheme.size()); 
	for (int i=0; i<(int)(block_scheme.size()); i++)
	{
		// B_matrix[i] = [e_{block_scheme[i](0)}, e_{block_scheme[i](1)}, ...]
		// block_scheme[i]: TIndex containing the dimensions of the i-th block
		B_matrix[i]=Zeros(x.data.dim, block_scheme[i].size); 
		B_matrix[i].Insert(block_scheme[i], TIndex(0,block_scheme[i].size-1), Identity(block_scheme[i].size)); 	
	}

	// double accP=1.0- exp(log(1.0-0.25)/(double)block_scheme.size());
	double accP = 0.234;
	BlockAdaptive(x, B_matrix, accP, period, max_period, if_eejump); 

	if (block_file_name.empty())
		return true; 
	else 
	{
		bool return_value = WriteBlocks(block_file_name); 
		return return_value; 
	}
}

bool CMetropolis::AdaptiveAfterSimulation_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &sample_file_name, const string &block_file_name, bool if_eejump, const string &block_scheme_file) 
// 	variance matrix estimated from simulation with scaled columns / 1 block
// Here filename refers to the file where samples are stored to. 
// So we first load these samples, and then calculate mean and variance. 
// finally blocks will be written into a file so that simulation can start from there
{	
	vector<CSampleIDWeight> Y; 
	if (!LoadSampleFromFile(sample_file_name, Y) )
		return false; 		
	
	block_scheme = ReadBlockScheme(block_scheme_file); 
	if (block_scheme.empty())
		block_scheme.push_back(TIndex(0, Y[0].data.dim-1)); 

	vector<TDenseMatrix>B_matrix(block_scheme.size()); 
	TDenseVector mean, d_vector; 
	TDenseMatrix variance, U_matrix, V_matrix, D_matrix; 
	for (int i_block=0; i_block<(int)block_scheme.size(); i_block++)
	{
		// Compute variance and mean
		mean = Zeros(block_scheme[i_block].size);
		variance = Zeros(block_scheme[i_block].size, block_scheme[i_block].size); 
		for (int i=0; i<(int)Y.size(); i++)
		{
			mean = mean + Y[i].data.SubVector(block_scheme[i_block]);  
			variance = variance + Multiply(Y[i].data.SubVector(block_scheme[i_block]),Y[i].data.SubVector(block_scheme[i_block]));
		}
	
		mean = (1.0/(double)Y.size()) * mean; 
		variance = (1.0/(double)Y.size())*variance - Multiply(mean, mean); 
		variance = 0.5 * (variance+Transpose(variance)); 
 
		SVD(U_matrix, d_vector, V_matrix, variance); 
		for (int i=0; i<d_vector.dim; i++)
			d_vector[i] = sqrt(d_vector[i]); 
		D_matrix = DiagonalMatrix(d_vector); 
		U_matrix = U_matrix *D_matrix; 

		B_matrix[i_block]=Zeros(Y[0].data.dim, block_scheme[i_block].size); 
		B_matrix[i_block].Insert(block_scheme[i_block], TIndex(0,block_scheme[i_block].size-1), U_matrix); 
	}
	CSampleIDWeight x=adaptive_start_point; 

	// double accP=1.0- exp(log(1.0-0.25)/(double)block_scheme.size());
	double accP = 0.234; 
	BlockAdaptive(x, B_matrix, accP, period, max_period, if_eejump); 

	if (block_file_name.empty())
		return true; 
	else 
	{
		bool return_value = WriteBlocks(block_file_name); 
		return return_value; 
	}
}

// bool CMetropolis::AdaptiveAfterSimulation_WeightedSampling_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const vector<CSampleIDWeight> &Y, const vector<double> &weight, const string &block_file_name, bool if_eejump, const string &block_scheme_file)
bool CMetropolis::AdaptiveAfterSimulation_WeightedSampling_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, std::vector<TDenseMatrix> &B_matrix, const string &block_file_name, bool if_eejump)
{
	/*block_scheme = ReadBlockScheme(block_scheme_file); 
	if (block_scheme.empty())
		block_scheme.push_back(TIndex(0, Y[0].data.dim-1)); 

	vector<TDenseMatrix>B_matrix(block_scheme.size()); 
	TDenseVector mean, d_vector; 
	TDenseMatrix variance, U_matrix, V_matrix, D_matrix; 
	for (int i_block=0; i_block<(int)block_scheme.size(); i_block++)
	{
		// Compute variance and mean
		mean = Zeros(block_scheme[i_block].size);
		variance = Zeros(block_scheme[i_block].size, block_scheme[i_block].size); 
		for (int i=0; i<(int)Y.size(); i++)
		{
			mean = mean + weight[i] * Y[i].data.SubVector(block_scheme[i_block]);  
			variance = variance + weight[i] * Multiply(Y[i].data.SubVector(block_scheme[i_block]),Y[i].data.SubVector(block_scheme[i_block]));
		}
	
		variance = variance - Multiply(mean, mean); 
		variance = 0.5 * (variance+Transpose(variance)); 
 
		SVD(U_matrix, d_vector, V_matrix, variance); 
		for (int i=0; i<d_vector.dim; i++)
			d_vector[i] = sqrt(d_vector[i]); 
		D_matrix = DiagonalMatrix(d_vector); 
		U_matrix = U_matrix *D_matrix; 

		B_matrix[i_block]=Zeros(Y[0].data.dim, block_scheme[i_block].size); 
		B_matrix[i_block].Insert(block_scheme[i_block], TIndex(0,block_scheme[i_block].size-1), U_matrix); 
	}*/
	CSampleIDWeight x=adaptive_start_point; 

	// double accP=1.0- exp(log(1.0-0.25)/(double)block_scheme.size());
	double accP = 0.234; 
	BlockAdaptive(x, B_matrix, accP, period, max_period, if_eejump); 

	if (block_file_name.empty())
		return true; 
	else 
	{
		bool return_value = WriteBlocks(block_file_name); 
		return return_value; 
	}
}

vector<TIndex> ReadBlockScheme(const string &file_name)
{
	// Block scheme file is in binary format
	// number of blocks
	// number of indices, indices
	ifstream input; 
	input.open(file_name.c_str(), ios::in|ios::binary); 
	if (!input.is_open())
		return vector<TIndex>(0);  
	// number of blocks
	int n_blocks; 
	input.read((char*)&(n_blocks), sizeof(int)); 
	vector<TIndex>block_scheme(n_blocks);
	
	for (int i_block=0; i_block<n_blocks; i_block++)
	{
		// number of indices, and then indices
		int n_indices; 
		input.read((char*)&(n_indices), sizeof(int)); 
		int index; 
		for (int i_index=0; i_index<n_indices; i_index++)
		{
			input.read((char*)&(index), sizeof(int)); 
			block_scheme[i_block] += index; 
		}
	} 
	input.close(); 
	return block_scheme; 
}

std::vector<TDenseMatrix> GetBlockMatrix_WeightedSampling(const std::vector<CSampleIDWeight> &Y, const std::vector<double> &weight, const std::vector<TIndex> &block_scheme)
{
	std::vector<TDenseMatrix>B_matrix(block_scheme.size()); 
        TDenseVector mean, d_vector;
        TDenseMatrix variance, U_matrix, V_matrix, D_matrix;
        for (int i_block=0; i_block<(int)block_scheme.size(); i_block++)
        {
                mean = Zeros(block_scheme[i_block].size);
                variance = Zeros(block_scheme[i_block].size, block_scheme[i_block].size);
                for (int i=0; i<(int)Y.size(); i++)
                {
                        mean = mean + weight[i] * Y[i].data.SubVector(block_scheme[i_block]);
                        variance = variance + weight[i] * Multiply(Y[i].data.SubVector(block_scheme[i_block]),Y[i].data.SubVector(block_scheme[i_block]));
                 }

                variance = variance - Multiply(mean, mean);
                variance = 0.5 * (variance+Transpose(variance));

                SVD(U_matrix, d_vector, V_matrix, variance);
                for (int i=0; i<d_vector.dim; i++)
                        d_vector[i] = sqrt(d_vector[i]);
                D_matrix = DiagonalMatrix(d_vector);
                U_matrix = U_matrix *D_matrix;
                B_matrix[i_block]=Zeros(Y[0].data.dim, block_scheme[i_block].size);
                B_matrix[i_block].Insert(block_scheme[i_block], TIndex(0,block_scheme[i_block].size-1), U_matrix);
       }
        return B_matrix;
}

bool WriteBMatrixFile(const string &file_name, std::vector<TDenseMatrix> &B_matrix)
{
	fstream output; 
	output.open(file_name.c_str(), ios::out|ios::binary);
        if (!output.is_open())
		return false; 
	int n_blocks = (int)(B_matrix.size()); 
	output.write((char *)(&n_blocks), sizeof(int)); // number of blocks
	for (int i=0; i<n_blocks; i++)
		B_matrix[i].WriteBinary(output); 
	output.close(); 
	return true; 
}

std::vector<TDenseMatrix> ReadBMatrixFile(const string &file_name)
{
	fstream input; 
	input.open(file_name.c_str(), ios::in|ios::binary); 
	if (!input.is_open())
		return std::vector<TDenseMatrix>(0); 

	int n_blocks; 
	input.read((char*)&(n_blocks), sizeof(int)); 
	std::vector<TDenseMatrix> BMatrix(n_blocks); 
	for (int i=0; i<n_blocks; i++)
		BMatrix[i].ReadBinary(input); 
	input.close(); 
	return BMatrix; 
}
