#include "CMetropolis.h"
#include <cmath>
#include "dw_rand.h"
#include <fstream>

bool CMetropolis::RandomBlockRandomWalkMetropolis(double &log_posterior_y, CSampleIDWeight &y, const CSampleIDWeight &initial_x)
// Given
// 	inital_x:	initial value, vector
// Results:
// 	y:	new vector according to random walk-metropolis using random blocks of directios and scales
// 	log_posterior_y:	log_posterior(y)
// Returns:
// 	true:	if y is a new sample (different than initial_x)
// 	false:	if y is the same as initial_x
// Involving:
// 	random_block_assignments: random_block_assignments[i] contains the dimensions that are assigned into the i-th random block
// 	random_blocks:	random_blocks[i] contains the direction of the i-th dimension when it is in a block by itself
// 	random_block_scales:	random_block_scale[i][j]: scale factor when the i-th dimension is in a block of size (j+1)
{
	size_t n_blocks = random_block_assignments.size(); 	// number of blocks
	
	CSampleIDWeight x=initial_x, increment; 
	double log_previous = model->log_posterior_function(x), log_current; 	
	bool if_new_sample = false; 
	for (unsigned int i=0; i<n_blocks; i++)
	{
		size_t block_size = random_block_assignments[i].size();	// size of the current block
		increment.data.Zeros(x.data.dim); 
		for (unsigned int j=0; j<block_size; j++)
		{
			unsigned int d_index = random_block_assignments[i][j]; 
			increment.data = increment.data + random_blocks[d_index] * random_block_scales[d_index][block_size-1] * dw_gaussian_rnd();  
		}
		y.data = x.data + increment.data; 
		log_current = model->log_posterior_function(y); 
		if (log_current - log_previous >= log(dw_uniform_rnd() ) )
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

void CMetropolis::RandomBlockAdaptive(const CSampleIDWeight &adaptive_start_point, double target_ratio, size_t period, size_t max_period)
{
	size_t n_blocks = random_block_assignments.size(); 
	unsigned int n_draws = 0; 
	vector<unsigned int> n_accepted(n_blocks,0); 
	vector<unsigned int> begin_draws(n_blocks,0); 
	vector<unsigned int> begin_jumps(n_blocks,0); 
	vector<unsigned int> periods(n_blocks,period); 
	vector<unsigned int> end_draws(n_blocks,period); 

	TDenseVector previous_ratio(n_blocks,0.0); 
	TDenseVector scale(n_blocks,1.0); 
	TDenseVector best_scale(n_blocks,1.0); 
	TDenseVector low_scale(n_blocks,-1.0); 
	TDenseVector low_jump_ratio(n_blocks,-1.0); 
	TDenseVector high_scale(n_blocks,-1.0); 
	TDenseVector high_jump_ratio(n_blocks,-1.0); 	

	double log_mid = log(target_ratio);
	double lower_bound = exp(log_mid/0.2);
	double upper_bound = exp(log_mid/5.0); 

	CSampleIDWeight x=adaptive_start_point, y, increment; 
	double log_previous = model->log_posterior_function(x), log_current, new_scale, diff; 
	bool done = false; 
	unsigned int check = period; 
	while (!done)	
	{
		for (unsigned int i=0; i<n_blocks; i++)
		{
			// draw metropolis
			increment.data.Zeros(x.data.dim); 
			for (unsigned int j=0; j<random_block_assignments[i].size(); j++)
			{
				unsigned int d_index = random_block_assignments[i][j]; 
				increment.data = increment.data + random_blocks[d_index]*scale[i]*dw_gaussian_rnd(); 
			}
			y.data = x.data+increment.data; 
			log_current = model->log_posterior_function(y); 
			if (log_current-log_previous >= log(dw_uniform_rnd() ) )
			{
				x = y; 
				log_previous = log_current; 
				n_accepted[i] ++; 
			}
		}
		n_draws ++; 

		// rescale 
		if (n_draws == check )
		{
			done = true; 
			for (unsigned int i=0; i<n_blocks; i++)
			{
				// jump ratio
				previous_ratio[i] = (double)(n_accepted[i]-begin_jumps[i])/(double)(n_draws-begin_draws[i]); 
				// recompute scale?
				if (end_draws[i] <= n_draws)
				{
					// set new low or high bounds
					if (previous_ratio[i] < target_ratio)
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
							new_scale = (log_mid/log(previous_ratio[i]))*high_scale[i]; 
					}
					else if (high_jump_ratio[i] < 0.0)
					{
						best_scale[i] = scale[i]; 
						if (previous_ratio[i] < lower_bound)
							new_scale = 0.2*low_scale[i]; 
						else 
							new_scale = (log_mid/log(previous_ratio[i]))*low_scale[i]; 
					}
					else 
					{
						diff = high_jump_ratio[i]-low_jump_ratio[i]; 
						if (diff > 1.0e-6)
							new_scale = ((target_ratio-low_jump_ratio[i])*low_scale[i] + (high_jump_ratio[i]-target_ratio)*high_scale[i])/diff;
						else 
							new_scale = 0.5 * (low_scale[i]+high_scale[i]); 
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

	random_block_scales.resize(random_blocks.size()); 
	for (unsigned int i=0; i<n_blocks; i++)
	{
		size_t block_size = random_block_assignments[i].size(); 
		for (unsigned int j=0; j<block_size; j++)
		{
			unsigned int d_index = random_block_assignments[i][j]; 
			random_block_scales[d_index][block_size-1] = scale[i]; 
		}
	}
}

void CMetropolis::FourPassRandomBlockAdaptive(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, size_t avg_block_size)
// first pass: 
// 	identity variance matrix/n blocks
// second pass
// 	diagonal variance matrix/1 block
// third pass 
// 	simulate from second pass block structure 
// 	variance matrix from simulation / n blocks
// 	fill random_blocks
// fourth pass 
// 	random_block_assignment
// 	RandomBlockAdaptive
{
	CSampleIDWeight x=adaptive_start_point; 

	// Copy blocks to blocks-backup, because blocks will be used later	
	vector<TDenseMatrix> blocks_backup(blocks.size()); 
	for (unsigned int i=0; i<blocks.size(); i++)
		blocks_backup[i].CopyContent(blocks[i]); 

	// first pass
	size_t n_blocks = x.data.dim; 
	vector<TDenseMatrix> B_matrix(n_blocks); 
	TDenseMatrix I_matrix = Identity(x.data.dim); 
	for (unsigned int i=0; i<n_blocks; i++)
		B_matrix[i] = ColumnMatrix(ColumnVector(I_matrix,i)); 
	BlockAdaptive(x, B_matrix, 0.25, period, max_period);

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
	// burn in
	double log_posterior_y;
        CSampleIDWeight y;
        for (unsigned int t=0; t<burn_in; t++)
        {
                if (BlockRandomWalkMetropolis(log_posterior_y, y, x) )
                        x=y; 
        }

        vector<CSampleIDWeight> Y_simulation(n_draws);
	// simulation
	for (unsigned int i=0; i<n_draws; i++)
        {
                for (unsigned int t=0; t<thin; t++)
                {
                        if (BlockRandomWalkMetropolis(log_posterior_y, y, x) )
                                x=y; 
                }
                Y_simulation[i]=y;
        }
	
	// calculate variance and mean
	TDenseVector mean(y.data.dim, 0.0); 
	TDenseMatrix variance(y.data.dim,y.data.dim,0.0), U_matrix, V_matrix, D_matrix;
        TDenseVector d_vector;
        for (unsigned int i=0; i<n_draws; i++)
	{
		mean = mean + Y_simulation[i].data; 
                variance = variance + Multiply(Y_simulation[i].data,Y_simulation[i].data);
	}
	
	mean = (1.0/(double)n_draws)*mean; 
        variance = (0.5/(double)n_draws)*(variance+Transpose(variance)); 
	variance = variance - Multiply(mean,mean); 

        SVD(U_matrix, d_vector, V_matrix, variance);
        D_matrix = DiagonalMatrix(d_vector);
        U_matrix = U_matrix *D_matrix;

	// third pass: n blocks
	n_blocks = x.data.dim;
        B_matrix.resize(n_blocks);
        for (unsigned int i=0; i<n_blocks; i++)
                B_matrix[i] = ColumnMatrix(ColumnVector(U_matrix,i));
        BlockAdaptive(x, B_matrix, 0.25, period, max_period);
	// blocks should contain n_blocks matrices, each of them being a x.dim-by-1 matrix
	
	// forth pass: random block
	random_blocks.resize(x.data.dim); 
	random_block_scales.resize(x.data.dim); 
	for (unsigned int i=0; i<x.data.dim; i++)
	{
		random_blocks[i] = ColumnVector(blocks[i],0); 
		random_block_scales[i] = TDenseVector(x.data.dim,1.0); 
	}
	AssignDimensionsToRandomBlocks(x.data.dim, avg_block_size); 
	RandomBlockAdaptive(x, 0.25, period, max_period);

	// Copy back blocks
	blocks.resize(blocks_backup.size()); 
	for (unsigned int i=0; i<blocks_backup.size(); i++) 
		blocks[i].CopyContent(blocks_backup[i]); 
}

bool CMetropolis::OnePassRandomBlockAdaptive(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t avg_block_size, const string &file_name)
// Assume random_blocks (directions) and random_block_scales (scales) have already been obtained at the third pass of FourPassAdaptive and are stored in a file
// So we only need to load these parameters in, and then do RandomBlockAdaptive in one pass
// In practice, we can try CMetropolis::OnePassRandomBlockAdaptive first. If it fails (because the file does not exist or file loading fails), then we call CMetropolis::FourPassRandomBlockAdaptive
{
	if (!Read_RandomBlocks_RandomBlockScales(file_name))
	{
		cerr << "CMetropolis::OnePassRandomBlockAdaptive : error in loading " << file_name << endl; 
		return false; 
	}
	AssignDimensionsToRandomBlocks(adaptive_start_point.data.dim, avg_block_size);
        RandomBlockAdaptive(adaptive_start_point, 0.25, period, max_period);
	return true; 
}

void CMetropolis::AssignDimensionsToRandomBlocks(size_t n, size_t avg_block_size)
{
	unsigned int *index_array = new unsigned int[n]; 
	for (unsigned int i=0; i<n; i++)
		index_array[i] = i; 
	dw_random_shuffle(index_array, n, sizeof(unsigned int)); 

	size_t n_blocks = (size_t)ceil((double)n/(double)avg_block_size); 
	random_block_assignments.resize(n_blocks); 
	unsigned int counter = 0; 
	for (unsigned int i=0; i<n_blocks; i++)
	{
		random_block_assignments[i].resize(avg_block_size < n-counter ? avg_block_size : n-counter); 
		for (unsigned int j=0; j<random_block_assignments[i].size(); j++)
		{
			random_block_assignments[i][j] = index_array[counter]; 
			counter ++; 
		}
	}
}

bool CMetropolis::Write_RandomBlocks_RandomBlockScales(const string &file_name)
{
	ofstream oFile(file_name.c_str(), ios::out|ios::binary); 
	if (!oFile)
		return false; 
	size_t n = random_blocks.size(); 
	oFile.write((char*)(&n),sizeof(size_t)); 
	for (unsigned int i=0; i<n; i++)
	{
		size_t m=random_blocks[i].dim; 
		oFile.write((char*)(&m),sizeof(size_t)); 
		for (unsigned int j=0; j<m; j++)
		{
			double data = random_blocks[i](j); 
			oFile.write((char*)(&data),sizeof(double)); 
		}
	}
	n = random_block_scales.size(); 
	oFile.write((char*)(&n), sizeof(size_t)); 
	for (unsigned int i=0; i<n; i++)
	{
		size_t m = random_block_scales[i].dim; 
		oFile.write((char*)(&m), sizeof(size_t)); 
		for (unsigned int j=0; j<m; j++)
		{
			double data = random_block_scales[i](j); 
			oFile.write((char*)(&data),sizeof(double)); 
		}
	}
	oFile.close(); 
	return true; 
}

bool CMetropolis::Read_RandomBlocks_RandomBlockScales(const string &file_name)
{
	ifstream iFile(file_name.c_str(), ios::binary|ios::in); 
	if (!iFile)
		return false; 
	size_t n, m; 
	double data; 
	iFile.read((char*)(&n), sizeof(size_t)); 
	random_blocks.resize(n); 
	for (unsigned int i=0; i<n; i++)
	{
		iFile.read((char*)(&m), sizeof(size_t)); 
		random_blocks[i].Resize(m); 
		for (unsigned int j=0; j<m; j++)
		{
			iFile.read((char*)(&data),sizeof(double)); 
			random_blocks[i](j) = data; 
		}
	}
	
        iFile.read((char*)(&n), sizeof(size_t));
	random_block_scales.resize(n); 
        for (unsigned int i=0; i<n; i++)
        {
		iFile.read((char*)(&m),sizeof(size_t)); 
		random_block_scales[i].Resize(m); 
                for (unsigned int j=0; j<m; j++)
                {
			iFile.read((char*)(&data),sizeof(double)); 
			random_block_scales[i](j) = data; 
                }
        }

	iFile.close(); 
	return true; 
}
