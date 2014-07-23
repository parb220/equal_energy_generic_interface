#include <cmath>
#include "prcsn.h"
#include "dw_rand.h"
#include "dw_dense_matrix.hpp"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CMetropolis.h"
#include "CEquiEnergyModel.h"

double CEquiEnergyModel::GMM_Simulation(int simulation_length, const std::vector<TDenseVector> &gmm_mean, const std::vector<TDenseMatrix> &gmm_covariance_sqrt)
{
	bool if_bounded_old = if_bounded;
        if_bounded = false;     // temperarily set if_bounded as false so hill-climbing is with respect to original model

	double log_posterior; 
	int N = (int)gmm_mean.size(), n; 
	int dim = current_sample.data.dim;
	int nAccpt = 0;  
	TDenseVector x(dim,0.0);  
	for (int i=0; i<simulation_length; i++)
	{
		for (int j=0; j<parameter->thin; j++)
		{
			n = dw_uniform_int(N); 
			x = gmm_covariance_sqrt[n]*x.RandomNormal()+gmm_mean[n];
			log_posterior = log_posterior_function(x.vector, x.dim); 
			if (log_posterior > MINUS_INFINITY)
			{
				CSampleIDWeight sample(x, (int)(time(NULL)-timer_when_started), log_posterior, true); 
 				SaveSampleToStorage(sample);
				nAccpt ++; 
			}
		} 
	}

	if_bounded = if_bounded_old;
	return (double)nAccpt/(double)(simulation_length); 
}

