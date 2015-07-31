#include <cmath>
#include <iomanip>
#include "prcsn.h"
#include "dw_rand.h"
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

void *MinusLogPosterior_NPSOL::function(int *mode, int *n, double *x, double *f, double *g, int *nstate)
{
	*f = -model->log_posterior_function(x,*n); 
}

double CEquiEnergyModel::HillClimb_NPSOL(int nSolution, int optimization_iteration, int perturbation_iteration, double perturbation_scale, double scale, const TDenseVector &start_point)
{
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
		bl[i] = MINUS_INFINITY; 
		bu[i] = PLUS_INFINITY; 
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

	int nAccpt = 0; 

	gmm_mean.clear(); 
	gmm_covariance_sqrt.clear();
	gmm_covariance_sqrt_log_determinant.clear(); 
	gmm_covariance_sqrt_inverse.clear(); 
 
	double log_posterior_before_perturbation, log_posterior_after_perturbation, log_posterior_optimal; 
	TDenseVector peak(n,0.0), perturbation(n,0.0), x_plus_perturbation(n,0.0);  
	TDenseMatrix hessian_cholesky(n,n,0.0); 

	try {
		for (int i=0; i<nSolution; i++)
		{
			if (start_point.dim)
				memcpy(x, start_point.vector, sizeof(double)*start_point.dim); 
			else 
			{
				log_posterior_before_perturbation = MINUS_INFINITY; 
				while (log_posterior_before_perturbation <= MINUS_INFINITY)
				{
					if(!DrawParametersFromPrior(x))
						throw dw_exception("CEquiEnergyModel::HillClimb_NPSOL : DrawSampleFromPrior() error occurred");
					log_posterior_before_perturbation = log_posterior_function(x, n);
				}
			}

			for (int ii=0; ii<optimization_iteration; ii++)
			{
				perturbation.RandomNormal(); 
				// Perturbation
				for (int jj=0; jj<perturbation_iteration; jj++)
				{
					memcpy(x_plus_perturbation.vector, x, sizeof(double)*n); 
					x_plus_perturbation += perturbation_scale * perturbation; 
					log_posterior_after_perturbation = log_posterior_function(x_plus_perturbation.vector, n);
					if (log_posterior_after_perturbation > log_posterior_before_perturbation)
					{
						memcpy(x, x_plus_perturbation.vector, sizeof(double)*n);
						log_posterior_before_perturbation = log_posterior_after_perturbation; 
					}	
					perturbation = perturbation*0.5; 	
				}
				
				// Optimization
				npoptn_((char*)DERIVATIVE_LEVEL.c_str(), DERIVATIVE_LEVEL.length());
                       	 	npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length());
                        	npoptn_((char*)COLD_START.c_str(), COLD_START.length());
                        	npoptn_((char*)HESSIAN.c_str(), HESSIAN.length()); 
                        	npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, MinusLogPosterior_NPSOL::function, &inform, &iter, istate, c, cJac, clambda, &f, g, R, x, iw, &leniw, w, &lenw);
				
				if (ii==0 || -f > log_posterior_optimal)
				{
					log_posterior_optimal = -f; 
					memcpy(peak.vector, x, sizeof(double)*n);
                        		memcpy(hessian_cholesky.matrix, R, sizeof(double)*n*n); 
                        		hessian_cholesky.column_major = true;
				}
			}

			// Accept the optimal solution
			if (log_posterior_optimal > MINUS_INFINITY)
			{
				nAccpt ++; 
				// Mean 
				TDenseVector copy(peak.dim); 
				copy.CopyContent(peak); 
				gmm_mean.push_back(copy);  
				// Hessian = inverse of covariance matrix
				TDenseMatrix hessian = Transpose(hessian_cholesky)*hessian_cholesky; 
				hessian = 0.5*(hessian+Transpose(hessian)); 
				TDenseMatrix eVectorHessian; 
				TDenseVector eValueHessian; 
				Eig(eValueHessian, eVectorHessian, hessian);  
				// hessian = eVectorHessian * DiagonalMatrix(eValueHessian) * eVectorHessian'
				// covariance_matrix = eVectorHessian * Inverse(DiagonalMatrix(eValueHessian)) *eVectorHessian';
			
				TDenseVector diffusedEValue(eValueHessian.dim,0.0), diffusedEValueInverse(eValueHessian.dim,0.0); 
				double covariance_sqrt_log_determinant = 0.0; 
				for (int i=0; i<eValueHessian.dim; i++)
				{
					diffusedEValue[i] = sqrt(scale/eValueHessian[i]); 
					diffusedEValueInverse[i] = sqrt(eValueHessian[i]/scale); 
					covariance_sqrt_log_determinant += -0.5*log(eValueHessian[i]/scale); 
				}

				gmm_covariance_sqrt.push_back(eVectorHessian*DiagonalMatrix(diffusedEValue)*Transpose(eVectorHessian)); 
				gmm_covariance_sqrt_inverse.push_back(eVectorHessian*DiagonalMatrix(diffusedEValueInverse)*Transpose(eVectorHessian)); 
				gmm_covariance_sqrt_log_determinant.push_back(covariance_sqrt_log_determinant); 
			}
		}
	}	
	catch( dw_exception& e )
	{
		cout << e.what() << endl; 
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

	return (double)nAccpt/(double)nSolution; 
}

bool CEquiEnergyModel::WriteGaussianMixtureModelParameters(const string &file_name) const
{
	int number = (int)gmm_mean.size();
        if (number != (int)gmm_covariance_sqrt.size())
                return false;
        fstream oFile(file_name.c_str(), ios::out | ios::binary);
        if (!oFile)
                return false;

        oFile.write((char *)(&number), sizeof(int));
        for (int i=0; i<number; i++)
        {
                oFile.write((char *)(&gmm_mean[i].dim), sizeof(int));
                oFile.write((char *)(gmm_mean[i].vector), sizeof(double)*gmm_mean[i].dim);
                oFile.write((char *)(&gmm_covariance_sqrt[i].rows), sizeof(int));
                oFile.write((char *)(&gmm_covariance_sqrt[i].cols), sizeof(int));
                oFile.write((char *)(&gmm_covariance_sqrt[i].column_major), sizeof(bool));
                oFile.write((char *)(gmm_covariance_sqrt[i].matrix), sizeof(double)*gmm_covariance_sqrt[i].rows*gmm_covariance_sqrt[i].cols);
		oFile.write((char *)&(gmm_covariance_sqrt_log_determinant[i]), sizeof(double)); 
		oFile.write((char *)&(gmm_covariance_sqrt_inverse[i].rows), sizeof(int)); 
		oFile.write((char *)&(gmm_covariance_sqrt_inverse[i].cols), sizeof(int)); 
		oFile.write((char *)&(gmm_covariance_sqrt_inverse[i].column_major), sizeof(bool)); 
		oFile.write((char *)(gmm_covariance_sqrt_inverse[i].matrix), sizeof(double)*gmm_covariance_sqrt_inverse[i].rows*gmm_covariance_sqrt_inverse[i].cols);  
        }
        oFile.close();

        return true;
}

bool CEquiEnergyModel::WriteGaussianMixtureModelMeanAscii(const string &file_name) 
{
	int number = (int)gmm_mean.size();
	ofstream oFile(file_name.c_str());
        if (!oFile)
                return false;
	double log_prior, log_likelihood; 
	for (int i=0; i<number; i++)
	{
		log_prior = log_prior_function(gmm_mean[i].vector, gmm_mean[i].dim); 
		log_likelihood = log_likelihood_function(gmm_mean[i].vector, gmm_mean[i].dim); 
		oFile << setprecision(20) << log_likelihood + log_prior << "\t" << log_likelihood; 
		for (int j=0; j<gmm_mean[i].dim; j++)
			oFile << setprecision(20) << "\t" << gmm_mean[i][j] ; 
		oFile << endl; 
	}
	oFile.close(); 
	return true; 
}

bool CEquiEnergyModel::ReadGaussianMixtureModelParameters(const string &file_name)
{
        fstream iFile(file_name.c_str(), ios::in | ios::binary);
        if (!iFile)
                return false;
        int number, dim, rows, cols;
        bool column_major;
        iFile.read((char *)&number, sizeof(int));
	if ((int)gmm_mean.size() != number)
		gmm_mean.resize(number); 
	if ((int)gmm_covariance_sqrt.size() != number)
		gmm_covariance_sqrt.resize(number); 
	if ((int)gmm_covariance_sqrt_log_determinant.size() != number)
		gmm_covariance_sqrt_log_determinant.resize(number); 
	if ((int)gmm_covariance_sqrt_inverse.size() != number)
		gmm_covariance_sqrt_inverse.resize(number); 

        for(int j=0; j<number; j++)
        {
                iFile.read((char *)&dim, sizeof(int));
		gmm_mean[j].Resize(dim); 
                iFile.read((char *)gmm_mean[j].vector,sizeof(double)*dim);

                iFile.read((char *)&rows, sizeof(int));
                iFile.read((char *)&cols, sizeof(int));
                iFile.read((char *)&column_major, sizeof(bool));
		gmm_covariance_sqrt[j].Resize(rows, cols); 
                iFile.read((char *)gmm_covariance_sqrt[j].matrix, sizeof(double)*rows*cols);
                gmm_covariance_sqrt[j].column_major = column_major;

		iFile.read((char *)&(gmm_covariance_sqrt_log_determinant[j]), sizeof(double)); 

		iFile.read((char *)&rows, sizeof(int));
		iFile.read((char *)&cols, sizeof(int));
		iFile.read((char *)&column_major, sizeof(bool));
		gmm_covariance_sqrt_inverse[j].Resize(rows, cols); 
		iFile.read((char *)gmm_covariance_sqrt_inverse[j].matrix, sizeof(double)*rows*cols); 
		gmm_covariance_sqrt_inverse[j].column_major = column_major; 
        }
        iFile.close();
        return true;
}

bool CEquiEnergyModel::AggregateGaussianMixtureModelParameters(const string &file_name) 
{
	fstream iFile(file_name.c_str(), ios::in | ios::binary);
        if (!iFile)
                return false;
        int number, dim, rows, cols;
        bool column_major;
        iFile.read((char *)&number, sizeof(int));

        for(int j=0; j<number; j++)
        {
                iFile.read((char *)&dim, sizeof(int));
                TDenseVector read_vector(dim,0.0);
                iFile.read((char *)read_vector.vector,sizeof(double)*dim);
		gmm_mean.push_back(read_vector); 

                iFile.read((char *)&rows, sizeof(int));
                iFile.read((char *)&cols, sizeof(int));
                iFile.read((char *)&column_major, sizeof(bool));
                TDenseMatrix read_matrix(rows, cols, 0.0);
                iFile.read((char *)read_matrix.matrix, sizeof(double)*rows*cols);
                read_matrix.column_major = column_major;
		gmm_covariance_sqrt.push_back(read_matrix);

		double read_number; 
		iFile.read((char *)&read_number, sizeof(double)); 
		gmm_covariance_sqrt_log_determinant.push_back(read_number); 

		iFile.read((char *)&rows, sizeof(int));
                iFile.read((char *)&cols, sizeof(int));
                iFile.read((char *)&column_major, sizeof(bool));
		TDenseMatrix read_matrix_inverse(rows, cols, 0.0);
                iFile.read((char *)read_matrix_inverse.matrix, sizeof(double)*rows*cols);
                read_matrix_inverse.column_major = column_major;
 		gmm_covariance_sqrt_inverse.push_back(read_matrix_inverse); 
        }
        iFile.close();
        return true;
}

void CEquiEnergyModel::ClearGaussianMixtureModelParameters()
{
	gmm_mean.clear(); 
	gmm_covariance_sqrt.clear(); 
	gmm_covariance_sqrt_log_determinant.clear(); 
	gmm_covariance_sqrt_inverse.clear(); 
}

void CEquiEnergyModel::KeepOptimalGaussianMixtureModelParameters()
{
	TDenseVector optmial_mean = gmm_mean[0]; 
	TDenseMatrix optimal_covariance_sqrt = gmm_covariance_sqrt[0]; 
	double optimal_covariance_sqrt_log_determinant = gmm_covariance_sqrt_log_determinant[0]; 
	TDenseMatrix optimal_covariance_sqrt_inverse = gmm_covariance_sqrt_inverse[0]; 

	double optimal_value = log_posterior_function(gmm_mean[0].vector, gmm_mean[0].dim); 
	for (int i=1; i<(int)gmm_mean.size(); i++)
	{
		double value = log_posterior_function(gmm_mean[i].vector, gmm_mean[i].dim); 
		if (value > optimal_value)
		{
			optimal_value = value; 
			optmial_mean = gmm_mean[i];
			optimal_covariance_sqrt = gmm_covariance_sqrt[i];
			optimal_covariance_sqrt_log_determinant = gmm_covariance_sqrt_log_determinant[i];
			optimal_covariance_sqrt_inverse = gmm_covariance_sqrt_inverse[i]; 
		}
	}

	// Keep the optimal one and remove all the others
	gmm_mean[0] = optmial_mean; 
	gmm_covariance_sqrt[0] = optimal_covariance_sqrt; 
	gmm_covariance_sqrt_log_determinant[0] = optimal_covariance_sqrt_log_determinant; 
	gmm_covariance_sqrt_inverse[0] = optimal_covariance_sqrt_inverse; 
	for (int i=(int)(gmm_mean.size()-1); i>0; i--)
	{
		gmm_mean.pop_back(); 
		gmm_covariance_sqrt.pop_back(); 
		gmm_covariance_sqrt_log_determinant.pop_back(); 
		gmm_covariance_sqrt_inverse.pop_back(); 
	}
}
