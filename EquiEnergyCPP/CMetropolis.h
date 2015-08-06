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
	vector<TDenseMatrix> blocks; 	// directions of the block. blocks[i] is n-by-bi,  where n is the dimension of the sample, bi is the size of the i-th block. 
	vector<double> scale; 

public:

	// draw one sample
	bool BlockRandomWalkMetropolis(double &, CSampleIDWeight &, const CSampleIDWeight &x, int thin=1); 

	// IO: blocks
	bool WriteBlocks(const string &file_name); 
	bool ReadBlocks(const string &file_name);
	bool ReadBlockScheme(const string &file_name); 

	void GetBlockMatrix_WeightedSampling(const std::vector<CSampleIDWeight> &Y, const std::vector<double> &weight); 
	void SetScale(const vector<double> &_scale); 
	void SetScale(double _scale); 
	void SetBlocks(const vector<TDenseMatrix> &_blocks); 
	void SetBlockScheme(const vector<TIndex> &_block_scheme); 

public: 
	// Constructions destructions
	CMetropolis(CEquiEnergyModel* _model=NULL) : model(_model), block_scheme(), blocks(), scale() {} 
	~CMetropolis() {}
}; 

#endif
