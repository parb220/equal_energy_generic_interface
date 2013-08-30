#ifndef _CLASS_METHOPOLIS_
#define _CLASS_METHOPOLIS_

#include <vector>
#include "dw_dense_matrix.hpp"

/*extern "C" {
	#include "dw_switch.h"
}*/

using namespace std; 
class CSampleIDWeight; 
class CEquiEnergyModel;  

bool LoadSampleFromFile(const string &file_name, vector<CSampleIDWeight> &Y); 


class CMetropolis
{
public: 
	CEquiEnergyModel* model;	// pointers to CEquiEnergyModel
protected:
	vector<TDenseMatrix> blocks; 	// each block is a n-by-bi matrix, consisting of the directions and scales of the i-th block, where n is the dimension of the sample, bi is the size of the i-th block. 
	vector<vector<int> > random_block_assignments;	// random_block_assignments[i] consists of the dimensions that are assigned into the i-th random block. Each dimension is assigned to one and only one block
	vector<TDenseVector> random_blocks;	// random_blocks[i]:  the direction when the i-th dimension is in a block by itself
	vector<TDenseVector> random_block_scales;	// random_block_scales[i][j]: scale factor when the i-th dimension is in a block of size (j+1)
	void AssignDimensionsToRandomBlocks(size_t n, size_t avg_block_size); 

public:
	
	// learning blocks
	void BlockAdaptive(const CSampleIDWeight &adaptive_start_point, const vector<TDenseMatrix > &B, double target_ratio, size_t period, size_t max_period); 
	bool FourPassAdaptive_StartWithoutSampleFile(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, const string &block_file_name=string()); 
	bool FourPassAdaptive_StartWithSampleFile(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, const string &sample_file_name, const string &block_file_name=string()); 
	bool AdaptiveBeforeSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &block_file_name=string());
	bool AdaptiveAfterSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &sample_file_name, const string &block_file_name=string());

	// draw one sample
	bool BlockRandomWalkMetropolis(double &, CSampleIDWeight &, const CSampleIDWeight &x, size_t thin=1); 

	// Random block learning
	void RandomBlockAdaptive(const CSampleIDWeight &adaptive_start_point, double target_ratio, size_t period, size_t max_period); 
	bool FourPassRandomBlockAdaptive_StartWithoutSampleFile(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, size_t avg_block_size, const string &block_file_name=string()); 
	bool RandomBlockAdaptiveAfterSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t avg_block_size, const string &sample_file_name, const string &block_file_name=string()); 
	// Draw one sample using random blocks
	bool RandomBlockRandomWalkMetropolis(double &, CSampleIDWeight &, const CSampleIDWeight &x, size_t thin=1); 

	// IO: blocks, random_blocks, random_block_scales
	bool WriteBlocks(const string &file_name); 
	bool ReadBlocks(const string &file_name); 
	bool Write_RandomBlocks_RandomBlockScales(const string &file_name); 
	bool Read_RandomBlocks_RandomBlockScales(const string &file_name); 

public: 
	// Constructions destructions
	CMetropolis(CEquiEnergyModel* _model=NULL) : model(_model) {} 
	~CMetropolis() {}
	
}; 
#endif
