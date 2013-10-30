#include "stdafx.h"
#include "ap.h"
#include "dataanalysis.h"
#include "CEquiEnergyModel.h"

bool CEquiEnergyModel::k_means_clustering(size_t K, int level_index, vector<CSampleIDWeight> & centers) const
{
	vector <CSampleIDWeight> samples; 
	if (!storage->DrawAllSample(level_index, samples) || samples.size() < K)
		return false; 
	if (samples.size() == K)
	{
		if (centers.size() != K)
			centers.resize(K); 
		for (int ii=0; ii<(int)(centers.size()); ii++)
			centers[ii] = samples[ii]; 
		return true; 	
	}
	// raw data
	int n_rows = (int)samples.size(), n_cols = samples[0].data.dim; 
	double *raw_data = new double(n_rows*n_cols); 
	for (int ii=0; ii<n_rows; ii++)
		memcpy(raw_data+ii*n_cols, samples[ii].data.vector, n_cols*sizeof(double));
	samples.clear(); 
	
	// raw data converted to real_2d_array
	alglib::real_2d_array xy; 
	xy.setcontent(n_rows, n_cols, raw_data);  

	// k-means clustering
	alglib::clusterizerstate s;
    	alglib::kmeansreport rep;
	
	alglib::clusterizercreate(s);
	alglib::clusterizersetpoints(s, xy, n_rows, n_cols, 2);
	alglib::clusterizersetkmeanslimits(s, 5, 0);
	alglib::clusterizerrunkmeans(s, K, rep); 
	delete [] raw_data; 

	if (rep.terminationtype == 1)
	{
		for (int ii=0; ii<n_rows; ii++)
		{
			centers[ii].data.Resize(n_cols); 
			memcpy(centers[ii].data.vector, rep.c[ii], sizeof(double)*n_cols); 
		}
		return true; 
	}
	return false; 
}
