#include <vector>
#include <fstream>
#include "dw_dense_matrix.hpp"

using namespace std; 

bool WriteGaussianMixtureModelParameters(const string &file_name, const vector<TDenseVector> &gm_mean, const vector<TDenseMatrix> &gm_covariance_sqrt) 
{
	int number = (int)gm_mean.size(); 
	if (number != (int)gm_covariance_sqrt.size())
		return false; 
	fstream oFile(file_name.c_str(), ios::out | ios::binary); 
	if (!oFile)
		return false; 

	oFile.write((char *)(&number), sizeof(int));
        for (int i=0; i<number; i++)
        {
        	oFile.write((char *)(&gm_mean[i].dim), sizeof(int));
        	oFile.write((char *)(gm_mean[i].vector), sizeof(double)*gm_mean[i].dim);
        	oFile.write((char *)(&gm_covariance_sqrt[i].rows), sizeof(int));
        	oFile.write((char *)(&gm_covariance_sqrt[i].cols), sizeof(int));
		oFile.write((char *)(&gm_covariance_sqrt[i].column_major), sizeof(bool)); 
        	oFile.write((char *)(gm_covariance_sqrt[i].matrix), sizeof(double)*gm_covariance_sqrt[i].rows*gm_covariance_sqrt[i].cols);
        }
        oFile.close();

	return true; 
}

bool ReadGaussianMixtureModelParameters(const string &file_name, vector<TDenseVector> &gm_mean, vector<TDenseMatrix >&gm_covariance_sqrt)
{
	gm_mean.clear(); 
	gm_covariance_sqrt.clear(); 
	fstream iFile(file_name.c_str(), ios::in | ios::binary);
        if (!iFile)
		return false; 
	int number, dim, rows, cols; 
	bool column_major; 
	iFile.read((char *)&number, sizeof(int));
        for(int j=0; j<number; j++)
        {
        	iFile.read((char *)&dim, sizeof(int));
        	TDenseVector peak(dim,0.0);
        	iFile.read((char *)peak.vector,sizeof(int)*dim);
        	
		iFile.read((char *)&rows, sizeof(int));
        	iFile.read((char *)&cols, sizeof(int));
		iFile.read((char *)&column_major, sizeof(bool));
        	TDenseMatrix covariance_sqrt(rows,cols,0.0);
        	iFile.read((char *)covariance_sqrt.matrix, sizeof(int)*rows*cols);
		covariance_sqrt.column_major = column_major; 
        	
		gm_mean.push_back(peak);
        	gm_covariance_sqrt.push_back(covariance_sqrt);
        }
        iFile.close();
	return true; 

}
