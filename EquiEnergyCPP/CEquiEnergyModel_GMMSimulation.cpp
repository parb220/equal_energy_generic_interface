#include <cmath>
#include <fstream>
#include <vector>
#include "prcsn.h"
#include "dw_rand.h"
#include "dw_math.h"
#include "dw_dense_matrix.hpp"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CMetropolis.h"
#include "CEquiEnergyModel.h"

double CEquiEnergyModel::GMM_Simulation(int simulation_length)
{
	bool if_bounded_old = if_bounded;
        if_bounded = false;     // temperarily set if_bounded as false so hill-climbing is with respect to original model

	int nAccpt = 0;  
	for (int i=0; i<simulation_length; i++)
	{
		CSampleIDWeight sample; 
		if (StudentT_DrawSample(sample))
		// if (Cauchy_DrawSample(sample))
		// if (GMM_DrawSample(sample))
		{	
 			SaveSampleToStorage(sample);
			nAccpt ++; 
		} 
	}

	if_bounded = if_bounded_old;
	return (double)nAccpt/(double)(simulation_length*parameter->thin); 
}

bool CEquiEnergyModel::StudentT_DrawSample(CSampleIDWeight &y) 
{
	int dim = current_sample.data.dim;
	TDenseVector x(dim,0.0);
	for (int i=0; i<dim; i++)
		x[i] = dw_tdistribution_rnd(5.0); // degree of freedom by default = 5.0 
	x = gmm_covariance_sqrt[0]*x+gmm_mean[0];
	double log_posterior = log_posterior_function(x.vector, x.dim);
	if (log_posterior > MINUS_INFINITY)
	{
		y = CSampleIDWeight(x, timer_when_started, log_posterior, true);
		return true; 
	}
	else 
		return false; 
}

bool CEquiEnergyModel::Cauchy_DrawSample(CSampleIDWeight &y) 
{
	int dim = current_sample.data.dim;
	TDenseVector x(dim,0.0);
	for (int i=0; i<dim; i++)
		x[i] = dw_cauchy(1.0); 
	x = gmm_covariance_sqrt[0]*x+gmm_mean[0];
	double log_posterior = log_posterior_function(x.vector, x.dim);
	if (log_posterior > MINUS_INFINITY)
	{
		y = CSampleIDWeight(x, timer_when_started, log_posterior, true);
		return true; 
	}
	else 
		return false; 
}

bool CEquiEnergyModel::GMM_DrawSample(CSampleIDWeight &y) 
{
	int N = (int)gmm_mean.size(), n=dw_uniform_int(N);
	int dim = current_sample.data.dim;
	TDenseVector x(dim,0.0);
	x = gmm_covariance_sqrt[n]*x.RandomNormal()+gmm_mean[n];
	double log_posterior = log_posterior_function(x.vector, x.dim);
	if (log_posterior > MINUS_INFINITY)
	{
		y = CSampleIDWeight(x, timer_when_started, log_posterior, true);
		return true; 
	}
	else 
		return false; 
}

double CEquiEnergyModel::Cauchy_LogPDF(const CSampleIDWeight &x) const 
{
	if (gmm_mean.empty() || gmm_covariance_sqrt.size() != gmm_mean.size() || gmm_covariance_sqrt_log_determinant.size() != gmm_mean.size() || gmm_covariance_sqrt_inverse.size() != gmm_mean.size())
		return MINUS_INFINITY; 
	
	TDenseVector rotatedError = gmm_covariance_sqrt_inverse[0] * (x.data-gmm_mean[0]); 
	double log_element = gmm_covariance_sqrt_log_determinant[0] ; 
	for (int i=0; i<rotatedError.dim; i++)
		log_element += dw_log_cauchy(rotatedError[i], 1.0); 

	return log_element ; 
}

double CEquiEnergyModel::StudentT_LogPDF(const CSampleIDWeight &x) const 
{
	if (gmm_mean.empty() || gmm_covariance_sqrt.size() != gmm_mean.size() || gmm_covariance_sqrt_log_determinant.size() != gmm_mean.size() || gmm_covariance_sqrt_inverse.size() != gmm_mean.size())
		return MINUS_INFINITY; 
	
	TDenseVector rotatedError = gmm_covariance_sqrt_inverse[0] * (x.data-gmm_mean[0]); 
	double log_element = gmm_covariance_sqrt_log_determinant[0] ; 
	for (int i=0; i<rotatedError.dim; i++)
		log_element += log(dw_tdistribution_pdf(rotatedError[i], 5.0)); // by default, degree of freedom = 5.0 

	return log_element ; 
}

double CEquiEnergyModel::GMM_LogPDF(const CSampleIDWeight &x) const 
{
	if (gmm_mean.empty() || gmm_covariance_sqrt.size() != gmm_mean.size() || gmm_covariance_sqrt_log_determinant.size() != gmm_mean.size() || gmm_covariance_sqrt_inverse.size() != gmm_mean.size())
		return MINUS_INFINITY; 
	
	TDenseVector rotatedError = gmm_covariance_sqrt_inverse[0] * (x.data-gmm_mean[0]); 
	double log_element_1 = -x.data.dim*0.918938533204673-log((double)gmm_mean.size()) -1.0*gmm_covariance_sqrt_log_determinant[0]-0.5*InnerProduct(rotatedError, rotatedError); 
	for (int i=1; i<(int)(gmm_mean.size()); i++)
	{
		rotatedError = gmm_covariance_sqrt_inverse[i] * (x.data - gmm_mean[i]); 
		double log_element_2 = -x.data.dim*0.918938533204673-log((double)gmm_mean.size())-1.0*gmm_covariance_sqrt_log_determinant[i] - 0.5*InnerProduct(rotatedError, rotatedError); 
		log_element_1 = AddLogs(log_element_1, log_element_2); 	
	}
	return log_element_1 ; 
}

double CEquiEnergyModel::StudentT_LogRatio(const CSampleIDWeight &x, const CSampleIDWeight &y) const
{
	double log_pdf_x = StudentT_LogPDF(x); 
	double log_pdf_y = StudentT_LogPDF(y); 
	if (log_pdf_x > MINUS_INFINITY && log_pdf_y > MINUS_INFINITY)
		return log_pdf_x - log_pdf_y; 
	else 
		return MINUS_INFINITY; 
}

double CEquiEnergyModel::Cauchy_LogRatio(const CSampleIDWeight &x, const CSampleIDWeight &y) const
{
	double log_pdf_x = Cauchy_LogPDF(x); 
	double log_pdf_y = Cauchy_LogPDF(y); 
	if (log_pdf_x > MINUS_INFINITY && log_pdf_y > MINUS_INFINITY)
		return log_pdf_x - log_pdf_y; 
	else 
		return MINUS_INFINITY; 
}

double CEquiEnergyModel::GMM_LogRatio(const CSampleIDWeight &x, const CSampleIDWeight &y) const
{
	double log_pdf_x = GMM_LogPDF(x); 
	double log_pdf_y = GMM_LogPDF(y); 
	if (log_pdf_x > MINUS_INFINITY && log_pdf_y > MINUS_INFINITY)
		return log_pdf_x - log_pdf_y; 
	else 
		return MINUS_INFINITY; 
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
