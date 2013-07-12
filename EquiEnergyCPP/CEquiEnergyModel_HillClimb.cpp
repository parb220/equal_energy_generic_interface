#include <cmath>
#include "dw_rand.h"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"

extern "C" {
	void npoptn_(char *, int);
	void npsol_(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A, double *bl, double *bu, void *funcon(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate), void *funobj(int *mode, int *n, double *x, double *f, double *g, int *nstate), int *inform, int *iter, int *istate, double *c, double *cJac, double *clamda, double *f, double *g, double *R, double *x, int *iw, int *leniw, double *w, int *lenw);
}

CEquiEnergyModel* MinusLogPosterior_NPSOL::model; 
void *MinusLogPosterior_NPSOL::function(int *mode, int *n, double *x, double *f, double *g, int *nstate)
{
	*f = model->log_posterior_function(x,(size_t)*n); 
}


double CEquiEnergyModel::log_posterior_function(const double *x, size_t n)
{
        double *old_x = new double[n];
	ConvertThetaToFreeParameters(target_model, old_x);
        ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) );
	
	ConvertFreeParametersToTheta(target_model, (double*)x);
	ConvertFreeParametersToQ(target_model, (double*)x+NumberFreeParametersTheta(target_model) );
	double log_posterior = LogPosterior_StatesIntegratedOut(target_model);
	ConvertFreeParametersToTheta(target_model, old_x);
        ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) );
	delete [] old_x; 
	
	return if_bounded ? -(-log_posterior>h_bound ? -log_posterior:h_bound)/t_bound : log_posterior; 
}

double CEquiEnergyModel::log_likelihood_function(const double *x, size_t n)
{
	double *old_x = new double[n];
        ConvertThetaToFreeParameters(target_model, old_x);
        ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) );
        
        ConvertFreeParametersToTheta(target_model, (double *)x);
        ConvertFreeParametersToQ(target_model, (double *)x+NumberFreeParametersTheta(target_model) );
	double log_likelihood = LogLikelihood_StatesIntegratedOut(target_model);

	// post old_x back to target_model 
	ConvertFreeParametersToTheta(target_model, old_x);
	ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) );
	delete [] old_x; 
	return log_likelihood; 
}

// HillClimb is always one on the original model. Therefore, if_bounded will be set as false temperoraly
// so that all calculation can be perfomed on the original model. After HillClimb is finished, if_bounded
// will be set as its original value. 
// Samples generated during HillClimb will be saved into storage but always at the level of number_energy_level
// (the extra level)
void CEquiEnergyModel::HillClimb(size_t nSolution, CStorageHead &storage, const CEESParameter &parameter)
{
	energy_level = parameter.number_energy_level; 
	bool if_bounded_old = if_bounded; 
	if_bounded = false;	// temperarily set if_bounded as false so all calculation is done on original model
 
	double *x_old = new double[current_sample.data.dim]; 
	ConvertThetaToFreeParameters(target_model, x_old); 
	ConvertQToFreeParameters(target_model, x_old+NumberFreeParametersTheta(target_model)); 
	
	const string COLD_START = string("Cold Start");
        const string NO_PRINT_OUT = string("Major print level = 0");
        const string DERIVATIVE_LEVEL = string("Derivative level = 0");

	// npsol 
	MinusLogPosterior_NPSOL::model = this; 
	int n = current_sample.data.dim; // sample dimension
	// linear constraints
	int nclin = 0; 	// number of linear constraints
	int ldA = nclin > 1 ? nclin : 1; 	// number of rows of A 
	double *A = new double[ldA*n]; 
	
	// nonlinear constraints	
	int ncnln = 0; 	// number of nonlinear constraints	
	int ldJ = ncnln > 1 ? ncnln : 1; 	// number of rows of cJac; 
	double *cJac = new double[ldJ*n]; 
	
	// R provides the upper-triangular Cholesky factor of the initial approximation
	// of the Hessian of the Lagrangian
	int ldR = n; 	// number of rows of R;
	double *R = new double[ldR*n];   
	
	// 
	int nctotal = n + nclin + ncnln; 
	int *istate = new int[nctotal]; // array indicating whether the constaits are effective
	double *bl = new double[nctotal]; // lower bounds for samples, A and cJac
	double *bu = new double[nctotal]; // upper bounds for samples, A and cJac
	for (int i=0; i<n; i++)
	{
		bl[i] = 0.0; 
		bu[i] = 100.0; 
	}
	double *clambda = new double [nctotal]; // lagragin parameters of the constraints	
	for (int i=0; i<nctotal; i++)
		clambda[i] = 0; 

	// work space
	int leniw = 3*n + nclin + 2*ncnln; 
	int *iw = new int [leniw]; 	// integer work space
	int lenw; 
	if (nclin == 0 && ncnln == 0)
		lenw = 20*n; 
	else if (ncnln == 0)
		lenw = 2*n*n + 20*n + 11*nclin; 
	else 
		lenw = 2*n*n + n*nclin + 2*n*ncnln + 20*n + 11*nclin + 21*ncnln;
	double *w = new double[lenw]; 	// floating work space

	// initial estimate of the solution 
	double *x = new double[n]; 	
	for (unsigned int i=0; i<n; i++)
		x[i]=dw_gaussian_rnd();  

	// returning value of npsol
	int inform, iter;
	double *c = new double[1]; // because ncnln = 0
	double f; 	// value of objective f(x) at the final iterate
	double *g = new double[n]; // objective gradient	   

	npoptn_((char*)DERIVATIVE_LEVEL.c_str(), DERIVATIVE_LEVEL.length());
        npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length());

	// test if LogPosterior_StatesIntegratedOut works
	for (unsigned int i=0; i<nSolution; i++)
	{
        	npoptn_((char*)COLD_START.c_str(), COLD_START.length());
		npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, MinusLogPosterior_NPSOL::function, &inform, &iter, istate, c, cJac, clambda, &f, g, R, x, iw, &leniw, w, &lenw); 

		CSampleIDWeight sample; 
		sample.data.Resize(n); 
		for (unsigned int j=0; j<n; j++)
			sample.data[i] = x[i]; 
		sample.id = (unsigned int)(time(NULL)-timer_when_started);    
		log_posterior_function(sample); 

		int binIndex = parameter.BinIndex(sample.weight, energy_level);
                storage.DepositSample(binIndex, sample);
	}

	delete [] g; 
	delete [] c;
	delete [] w; 
	delete [] iw; 
	delete [] istate;
	delete [] clambda; 
	delete [] bu; 
	delete [] bl;  
	delete [] R; 
	delete [] cJac; 
	delete [] A; 

	storage.finalize(parameter.BinIndex_Start(energy_level), parameter.BinIndex_End(energy_level));  
	storage.ClearDepositDrawHistory(parameter.BinIndex_Start(energy_level), parameter.BinIndex_End(energy_level));

	// recover to the old seeting
	ConvertFreeParametersToTheta(target_model, x_old); 
	ConvertFreeParametersToQ(target_model, x_old+NumberFreeParametersTheta(target_model)); 
	delete [] x_old;
	if_bounded = if_bounded_old; 
}
