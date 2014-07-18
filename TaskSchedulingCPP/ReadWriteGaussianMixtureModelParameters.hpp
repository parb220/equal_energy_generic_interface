#ifndef _RW_GMM_HEADER_
#define _RW_GMM_HEADER_

#include <string>
#include <vector>
#include "dw_dense_matrix.hpp"

bool WriteGaussianMixtureModelParameters(const string &file_name, const vector<TDenseVector> &gm_mean, const vector<TDenseMatrix> &gm_covariance_sqrt); 

bool ReadGaussianMixtureModelParameters(const string &file_name, vector<TDenseVector> &gm_mean, vector<TDenseMatrix >&gm_covariance_sqrt); 

#endif
