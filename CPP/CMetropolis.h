#ifndef _CLASS_METHOPOLIS_
#define _CLASS_METHOPOLIS_

#include <vector>
#include "CEquiEnergyModel.h"
#include "dw_dense_matrix.hpp"

extern "C" {
	#include "dw_switch.h"
}

using namespace std; 

class CMetropolis
{
protected:
	vector<TDenseMatrix> blocks; 	// each block is a n-by-bi matrix, consisting of the directions and scales of the i-th block, where n is the dimension of the sample, bi is the size of the i-th block. 
	vector<vector<unsigned int> > random_block_assignments;	// random_block_assignments[i] consists of the dimensions that are assigned into the i-th random block. Each dimension is assigned to one and only one block
	vector<TDenseVector> random_blocks;	// random_blocks[i]:  the direction when the i-th dimension is in a block by itself
	vector<TDenseVector> random_block_scales;	// random_block_scales[i][j]: scale factor when the i-th dimension is in a block of size (j+1)

	void AssignDimensionsToRandomBlocks(size_t n, size_t avg_block_size); 

public:
	CEquiEnergyModel *model;	// pointers to CEquiEnergyModel 
	
	// learning blocks
	void BlockAdaptive(const TDenseVector &adaptive_start_point, const vector<TDenseMatrix > &B, double target_ratio, size_t period, size_t max_period); 
	void FourPassAdaptive(const TDenseVector &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, const string &rb_file=string()); 

	// draw one sample
	bool BlockRandomWalkMetropolis(double &, TDenseVector &, const TDenseVector &x); 

	// Random block learning
	void RandomBlockAdaptive(const TDenseVector &adaptive_start_point, double target_ratio, size_t period, size_t max_period); 
	void FourPassRandomBlockAdaptive(const TDenseVector &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, size_t avg_block_size); 
	void OnePassRandomBlockAdaptive(const TDenseVector &adaptive_start_point, size_t period, size_t max_period, size_t avg_block_size); 
	// Draw one sample using random blocks
	bool RandomBlockRandomWalkMetropolis(double &, TDenseVector &, const TDenseVector &x); 

	// IO: blocks, random_blocks, random_block_scales
	bool WriteBlocks(const string &file_name); 
	bool ReadBlocks(const string &file_name); 
	bool Write_RandomBlocks_RandomBlockScales(const string &file_name); 
	bool Read_RandomBlocks_RandomBlockScales(const string &file_name); 
}; 
#endif
