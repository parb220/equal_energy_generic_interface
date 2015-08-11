#ifndef _CLASS_METHOPOLIS_
#define _CLASS_METHOPOLIS_

#include <vector>
#include "dw_dense_matrix.hpp"

using namespace std; 
class CSampleIDWeight; 
class CEquiEnergyModel;  


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
	void SimpleBlockAdaptive(const CSampleIDWeight &adaptive_start_point, const vector<TDenseMatrix> &B, int period, double lower_bound, double upper_bound, bool if_eejump=false); 
	
	bool AdaptiveBeforeSimulation_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, const string &block_file_name=string(), bool if_eejump=false, const string &block_scheme_file_name=string());

	bool AdaptiveAfterSimulation_WeightedSampling_OnePass(const CSampleIDWeight &adaptive_start_point, size_t period, size_t max_period, std::vector<TDenseMatrix> &B_matrix, const string &block_file_name=string(), bool if_eejump=false); 

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

std::vector<TIndex>ReadBlockScheme(const string &file_name); 
std::vector<TDenseMatrix> ReadBMatrixFile(const string &file_name); 
bool WriteBMatrixFile(const string &file_name, std::vector<TDenseMatrix> &B_matrix);
std::vector<TDenseMatrix> GetBlockMatrix_WeightedSampling(const std::vector<CSampleIDWeight> &Y, const std::vector<double> &weight, const std::vector<TIndex> &block_scheme); 

#endif
