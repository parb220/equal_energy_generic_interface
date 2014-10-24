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
	vector<TIndex> block_scheme; 	// assignments of the dimensions into the blocks. block_scheme[i] contains the dimensions assgined to the i-th block, which could be continuous or discontinues (general case)
	vector<TDenseMatrix> blocks; 	// directions and scales of the block. blocks[i] is n-by-bi,  where n is the dimension of the sample, bi is the size of the i-th block. 

public:
	
	// learning blocks
	void BlockAdaptive(const CSampleIDWeight &adaptive_start_point, const vector<TDenseMatrix > &B, double target_ratio, size_t period, size_t max_period, bool if_eejump =false); 
	bool FourPassAdaptive_StartWithoutSampleFile(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, const string &block_file_name=string(), bool if_eejump=false); 
	bool FourPassAdaptive_StartWithSampleFile(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, size_t n_draws, size_t burn_in, size_t thin, const string &sample_file_name, const string &block_file_name=string(), bool if_eejump=false); 
	bool AdaptiveBeforeSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &block_file_name=string(), bool if_eejump=false);
	bool AdaptiveAfterSimulation(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &sample_file_name, const string &block_file_name=string(), bool if_eejump=false);
	bool AdaptiveBeforeSimulation_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &block_file_name=string(), bool if_eejump=false, const string &block_scheme_file_name=string());
	bool AdaptiveAfterSimulation_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &sample_file_name, const string &block_file_name=string(), bool if_eejump=false, const string &block_scheme_file_name=string());

	bool AdaptiveAfterSimulation_WeightedSampling_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const vector<CSampleIDWeight> &samples, const vector<double> &weights, const string &block_file_name=string(), bool if_eejump=false, const string &block_scheme_file_name=string()); 

	// draw one sample
	bool BlockRandomWalkMetropolis(double &, CSampleIDWeight &, const CSampleIDWeight &x, size_t thin=1); 

	// IO: blocks
	bool WriteBlocks(const string &file_name); 
	bool ReadBlocks(const string &file_name);
	bool AggregateBlocksAndRemoveFiles(const vector<string> &read_file, const string &write_file_name); 

public: 
	// Constructions destructions
	CMetropolis(CEquiEnergyModel* _model=NULL) : model(_model) {} 
	~CMetropolis() {}
}; 

vector<TIndex>ReadBlockScheme(const string &file_name); 
#endif
