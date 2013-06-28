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
public:
	CEquiEnergyModel *model;	// pointers to CEquiEnergyModel 
	vector<TDenseMatrix> blocks; 
	
	// learning blocks
	void BlockAdaptive(const TDenseVector &adaptive_start_point, const vector<TDenseMatrix > &B, double target_ratio, int period, int max_period); 
	// draw one sample
	bool BlockRandomWalkMetropolis(double &, TDenseVector &, const TDenseVector &x); 
}; 
#endif
