#include <cmath>
#include "prcsn.h"
#include "dw_rand.h"
#include "dw_dense_matrix.hpp"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CMetropolis.h"
#include "CEquiEnergyModel.h"
#include "dw_dense_matrix.hpp"

extern "C" {
	void npoptn_(char *, int);
	void npsol_(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A, double *bl, double *bu, void *funcon(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate), void *funobj(int *mode, int *n, double *x, double *f, double *g, int *nstate), int *inform, int *iter, int *istate, double *c, double *cJac, double *clamda, double *f, double *g, double *R, double *x, int *iw, int *leniw, double *w, int *lenw);
	// #include "dw_csminwel.h"
}

CEquiEnergyModel* MinusLogPosterior_NPSOL::model; 
// CEquiEnergyModel* MinusLogPosterior_CSMINWEL::model; 

void *MinusLogPosterior_NPSOL::function(int *mode, int *n, double *x, double *f, double *g, int *nstate)
{
	*f = -model->log_posterior_function(x,*n); 
}

/*double MinusLogPosterior_CSMINWEL::function(double *x, int n, double **args, int *dims)
{
	double return_value = -model->log_posterior_function(x, n); 
	return return_value; 
}*/

// HillClimb is always one on the original model. Therefore, if_bounded will be set as false temperoraly
// so that all calculation can be perfomed on the original model. After HillClimb is finished, if_bounded
// will be set as its original value. 
// Samples generated during HillClimb will be saved into storage but always at the level of number_energy_level
// (the extra level)
double CEquiEnergyModel::HillClimb_NPSOL(int nSolution, std::vector<TDenseVector> &gm_mean, std::vector<TDenseMatrix> &gm_covariance_sqrt)
{
	// energy_level = parameter->number_energy_level; 
	bool if_bounded_old = if_bounded; 
	if_bounded = false;	// temperarily set if_bounded as false so hill-climbing is with respect to original model
 
	const string COLD_START = string("Cold Start");
        const string NO_PRINT_OUT = string("Major print level = 0");
        const string DERIVATIVE_LEVEL = string("Derivative level = 0");
	const string HESSIAN = string("Hessian = Yes"); 

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

	// returning value of npsol
	int inform, iter;
	double *c = new double[1]; // because ncnln = 0
	double f; 	// value of objective f(x) at the final iterate
	double *g = new double[n]; // objective gradient	   

	double max_log_posterior = MINUS_INFINITY; 

	gm_mean.clear(); 
	gm_covariance_sqrt.clear(); 
	try {
		for (int i=0; i<nSolution; i++)
		{
			if (!DrawParametersFromPrior(x)) 
				throw dw_exception("CEquiEnergyModel::HillClimb_NPSOL : DrawSampleFromPrior() error occurred");
			npoptn_((char*)DERIVATIVE_LEVEL.c_str(), DERIVATIVE_LEVEL.length());
        		npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length());
        		npoptn_((char*)COLD_START.c_str(), COLD_START.length());
			npoptn_((char*)HESSIAN.c_str(), HESSIAN.length()); 
			npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, MinusLogPosterior_NPSOL::function, &inform, &iter, istate, c, cJac, clambda, &f, g, R, x, iw, &leniw, w, &lenw); 
			if (f < PLUS_INFINITY)
			{
				/*CSampleIDWeight sample; 
				sample.data.Resize(n); 
				for (int j=0; j<n; j++)
					sample.data[j] = x[j]; 
				sample.DataChanged(); 
				sample.id = (int)(time(NULL)-timer_when_started);    
				log_posterior_function(sample); 

        			SaveSampleToStorage(sample);  */       	
				max_log_posterior = -f > max_log_posterior ? -f : max_log_posterior; 

				// Mean 
				TDenseVector peak(n,0.0); 
				memcpy(peak.vector, x, sizeof(double)*n);
				gm_mean.push_back(peak);  
				// Hessian (inverse of covariance matrix)
				TDenseMatrix hessian_cholesky(n,n,0.0); 
				memcpy(hessian_cholesky.matrix, R, sizeof(double)*n*n); 
				hessian_cholesky.column_major = true; 
				TDenseMatrix hessian = Transpose(hessian_cholesky)*hessian_cholesky; 
				hessian = 0.5*(hessian+Transpose(hessian)); 
				TDenseMatrix eigenVector; 
				TDenseVector eigenValue; 
				Eig(eigenValue, eigenVector, hessian); 
				for (int i=0; i<eigenValue.dim; i++)
					eigenValue[i] = sqrt(3.0/eigenValue[i]); // diffuse covariance matrix and take root
				TDenseMatrix covariance_sqrt = eigenVector * DiagonalMatrix(eigenValue); 
				gm_covariance_sqrt.push_back(covariance_sqrt); 
			}
		}
	}	
	catch(...)
	{
		if_bounded = if_bounded_old; 
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
		throw; 
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

	// storage->finalize(energy_level); 
	// storage->ClearDepositDrawHistory(energy_level);

	if_bounded = if_bounded_old; 
	return max_log_posterior; 
}

/*double CEquiEnergy_TState::HillClimb_CSMINWEL(size_t nSolution)
{
	energy_level = parameter->number_energy_level;
        bool if_bounded_old = if_bounded;
        if_bounded = false; 

	MinusLogPosterior_CSMINWEL::model = this;	

	int n = current_sample.data.dim;
	double *H=new double[n*n], *g=new double[n], *x=new double[n], crit = 1.0e-3, fh; 
	int nit=50, itct, fcount; 
	int retcodeh; 
	const double IniHCsminwel=1.0e-5; 
	double max_log_posterior = MINUS_INFINITY; 
	
	for (int iSolution=0; iSolution < nSolution; iSolution ++)
	{
		InitializeParameter(x, n); 
		for (int i=0; i<n; i++)
			for (int j=0; j<n; j++)
				H[i*n+j] = IniHCsminwel; 
		for (int i=0; i<n; i++)
			g[i] = 0.0; 
	
		dw_csminwel(MinusLogPosterior_CSMINWEL::function, x, n, H, g, NULL, &fh, crit, &itct, nit, &fcount, &retcodeh, NULL, NULL);
		if (retcodeh == 0)
		{
			CSampleIDWeight sample;
                        sample.data.Resize(n);
                        for (int i=0; i<n; i++)
                                sample.data[i] = x[i];
                        sample.DataChanged();
                        sample.id = (int)(time(NULL)-timer_when_started);
                        log_posterior_function(sample);

			SaveSampleToStorage(sample); 
			max_log_posterior = sample.weight > max_log_posterior ? sample.weight : max_log_posterior; 
		}
	}

	delete []x; 
	delete []g; 
	delete []H; 
	
	storage->finalize(energy_level); 
        storage->ClearDepositDrawHistory(energy_level); 
	if_bounded = if_bounded_old; 
	return max_log_posterior; 
}*/
