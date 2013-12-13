#include <cmath>
#include <algorithm>
#include "stdafx.h"
#include "ap.h"
#include "dw_rand.h"
#include "dataanalysis.h"
#include "CEquiEnergyModel.h"

bool CEquiEnergyModel::Initialize_MostDistant_WithinPercentile(size_t K, int level_index, vector<CSampleIDWeight > &starters, double percentile) const
{
        if (starters.size() != K)
                starters.resize(K);
	int nSample=0; 
	for (int ii=0; ii<storage->Number_Bin(level_index); ii++)
		nSample += (int)storage->GetNumberRecrod(level_index, ii); 
	if (nSample < K)
		return false; 
	
	vector<CSampleIDWeight> samples; 
	int nDraw = ceil(percentile*nSample) < K ? K : ceil(percentile*nSample); 
	int bin = 0;
       	while (bin < storage->Number_Bin(level_index) && samples.size() < nDraw)
       	{
               	storage->Draw_K_MostWeightSample(nDraw, level_index, bin, samples);
               	bin ++;
       	}
	if (nDraw == K)
	{
		for (int ii=0; ii<K; ii++)
			starters[ii] = samples[ii]; 
		return true; 
	}
	// 1st starter, the best sample so far
	starters[0] = samples[0]; 
	// k-th starter, the most distant to starters(0..k-1) 
	double most_distance, distance; 
	TDenseVector e; 
	for (int ii=1; ii<K; ii++)
	{
		// examining all the samples
		for(int jj=0; jj<(int)(samples.size()); jj++)
		{
			distance = 0.0; 
			// calculate the distance to starters[0, ..., ii-1]
			for (int kk=0; kk<ii; kk++)
			{
				e = samples[jj].data-starters[kk].data; 
				distance += Norm(e)*Norm(e); 
			}
			// use the sample that is fartherest to all existing centers
			if (jj == 0 || most_distance < distance)
			{
				most_distance = distance; 
				starters[ii] = samples[jj]; 
			}
		}
	}
	samples.clear(); 
	return true; 
}

bool CEquiEnergyModel::Initialize_MostDistant_WithinPercentileBand(size_t K, int level_index, vector<CSampleIDWeight > &starters, double percentile) const
{
        if (starters.size() != K)
                starters.resize(K);
	int nSample=0; 
	for (int ii=0; ii<storage->Number_Bin(level_index); ii++)
		nSample += (int)storage->GetNumberRecrod(level_index, ii); 
	if (nSample < K)
		return false; 
	
	vector<CSampleIDWeight> samples; 
	int nDraw = ceil(percentile*nSample) < K ? K : ceil(percentile*nSample); 
	int bin = 0;
       	while (bin < storage->Number_Bin(level_index) && samples.size() < nDraw)
       	{
               	storage->Draw_K_MostWeightSample(nDraw, level_index, bin, samples);
               	bin ++;
       	}
	if (nDraw == K)
	{
		for (int ii=0; ii<K; ii++)
			starters[ii] = samples[ii]; 
		return true; 
	}
	// 1st starter, the best sample so far
	starters[0] = samples[0]; 
	// k-th starter, the most distant to starters(0..k-1) within the range [being_index, end_index)
	int begin_index=1, end_index; 
	double most_distance, distance; 
	TDenseVector e; 
	for (int ii=1; ii<K; ii++)
	{
		end_index=(int)samples.size()*ii/K; 
		// examining the samples within the range [being_index, end_index)
		for(int jj=begin_index; jj<end_index; jj++)
		{
			distance = 0; 
			// calculate the distance to starters[0, ..., ii-1]
			for (int kk=0; kk<ii; kk++)
			{
				e = samples[jj].data-starters[kk].data; 
				distance += Norm(e)*Norm(e); 
			}
			// use the sample that is fartherest to all existing centers
			if (jj == begin_index || most_distance < distance)
			{
				most_distance = distance; 
				starters[ii] = samples[jj]; 
			}
		}
		begin_index = end_index; 
	}
	samples.clear(); 
	return true; 
}


bool CEquiEnergyModel::Initialize_KMeansClustering(size_t K, int level_index, vector<CSampleIDWeight> & centers) const
{
	if (centers.size() != K)
		centers.resize(K); 
	vector <CSampleIDWeight> samples; 
	if (!storage->DrawAllSample(level_index, samples) || samples.size() < K)
		return false; 
	if (samples.size() == K)
	{
		for (int ii=0; ii<(int)(centers.size()); ii++)
			centers[ii] = samples[ii]; 
		return true; 	
	}
	// raw data
	int n_rows = (int)samples.size(), n_cols = samples[0].data.dim; 
	double *raw_data = new double[n_rows*n_cols]; 
	for (int ii=0; ii<n_rows; ii++)
		memcpy(raw_data+ii*n_cols, samples[ii].data.vector, n_cols*sizeof(double));
	
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
		for (int ii=0; ii<K; ii++)
		{
			centers[ii].data.Resize(n_cols); 
			memcpy(centers[ii].data.vector, rep.c[ii], sizeof(double)*n_cols); 
		}
		return true; 
	}
	samples.clear(); 
	return false; 
}

bool CEquiEnergyModel::InitializeFromFile(const string &file_name)
{
        CSampleIDWeight x;
        ifstream input_file;
        input_file.open(file_name.c_str(), ios::binary|ios::in);
        if (!input_file)
                return false;
        read(input_file, &(x));
        input_file.close();
        Take_Sample_Just_Drawn_From_Storage(x);
        current_sample.DataChanged();
        log_posterior_function(current_sample);
        return true;
}

// A sample is randomly taken from a pool (with size desired_pool_size) of samples where the pool is formed by samples with higher log-posteriors. Note that samples with higher log-posterior values are stored in smaller-indexed bins. So, if the desired pool size is 10 while the size of the first bin is 100, then only the first bin will be used. In contrast, if the desired pool size is 100 while the total number of samples in the first 3 bins is barely greater than 100, then the first 3 bins will be used. 
bool CEquiEnergyModel::Initialize(size_t desiredPoolSize, int level_index)
{
        size_t N=storage->Number_Bin(level_index);
        int nSample_Total=0;
        vector<int>nSample_Bin(N,0);
        for (int bin_index=0;  bin_index<N; bin_index++)
        {
                nSample_Bin[bin_index] = storage->GetNumberRecrod(level_index, bin_index);
                nSample_Total += nSample_Bin[bin_index];
        }

        int nLumSum = 0, bin_index=0;
        int random_index = dw_uniform_int(desiredPoolSize < nSample_Total ? desiredPoolSize : nSample_Total);
        while (bin_index<N && !(random_index >= nLumSum && random_index < nLumSum + nSample_Bin[bin_index]) )
        {
                nLumSum += nSample_Bin[bin_index];
                bin_index++;
        }
        if (random_index >= nLumSum && random_index < nLumSum + nSample_Bin[bin_index])
        {
                if(storage->DrawSample(level_index, bin_index, current_sample))
                {
                        current_sample.id = (int)(time(NULL)-timer_when_started);
			// Because all samples stored in storage have had their log-posterior calculated and stored 
			// together with the sample values, there is no need to recalculate log-posterior at this moment
		return true;
                }
        }
        return false;
}

bool CEquiEnergyModel::InitializeWithBestSample(int level_index)
{
        int bin_index=0;
        while (bin_index<storage->Number_Bin(level_index) && !(storage->DrawMostWeightSample(level_index, bin_index, current_sample) ) )
                bin_index ++;
        if (bin_index >= storage->Number_Bin(level_index))
                return false;
        current_sample.id = (int)(time(NULL)-timer_when_started);
        return true;
}

bool CEquiEnergyModel::InitializeWith_Kth_BestSample(size_t K, int level_index)
{
        vector<CSampleIDWeight> sample;
        int bin=0;
        while (bin < storage->Number_Bin(level_index) && sample.size() < K)
        {
                storage->Draw_K_MostWeightSample(K, level_index, bin, sample);
                bin ++;
        }
	Take_Sample_Just_Drawn_From_Storage(sample.back());
        return true;
}

bool CEquiEnergyModel::Initialize_RandomlyPickFrom_K_BestSample(size_t K, int level_index)
{
        vector<CSampleIDWeight> sample;
        int bin = 0;
        while (bin < storage->Number_Bin(level_index) && sample.size() < K)
        {
                storage->Draw_K_MostWeightSample(K, level_index, bin, sample);
                bin ++;
        }
        //if (sample.size() < K)
        //      return false; 
        std::random_shuffle(sample.begin(), sample.end()); 
        int index = dw_uniform_int(K < sample.size() ? K : sample.size());
        Take_Sample_Just_Drawn_From_Storage(sample[index]);
        return true;
}

