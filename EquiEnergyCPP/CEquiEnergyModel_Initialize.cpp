#include <cmath>
#include <algorithm>
#include "dw_rand.h"
#include "dw_math.h"
#include "CEquiEnergyModel.h"

vector<double> CEquiEnergyModel::Reweight(const vector<CSampleIDWeight> &samples, int current_stage, int previous_stage)  const
{
	vector<double> log_weight(samples.size(),0.0); 
	// log_weight = weight*(1.0/t[current_stage]-1.0/([previous_stage]))
	for (int i=0; i<(int)(samples.size()); i++)
	{
		log_weight[i] = samples[i].weight  * parameter->lambda[current_stage]; 
		log_weight[i] = log_weight[i] - samples[i].weight * parameter->lambda[previous_stage] ; 
	}
	return log_weight; 
}

std::vector<CSampleIDWeight> CEquiEnergyModel::Initialize_WeightedSampling(const std::vector<CSampleIDWeight> &samples, int K, int stage_index) const
{
	if (!samples.empty())
	{
		std::vector<CSampleIDWeight> starters(K); 
		// Cumulative sum of importance weights
		vector<double> log_weight = Reweight(samples, stage_index-1, stage_index); 
		vector<double> log_weight_sum(samples.size()), weight_sum(samples.size()); 
		log_weight_sum[0] = log_weight[0];  
		for (int i=1; i<(int)samples.size(); i++)
			log_weight_sum[i] = AddLogs(log_weight_sum[i-1], log_weight[i]); 
		for (int i=0; i<(int)samples.size(); i++) // Normalize
			weight_sum[i] = exp(log_weight_sum[i] -log_weight_sum.back()); 

		for (int i=0; i<K; i++)
		{
			double random_number = dw_uniform_rnd(); 
			int position = std::lower_bound(weight_sum.begin(), weight_sum.end(), random_number)-weight_sum.begin(); 
			starters[i] = samples[position]; 	
		} 
		return starters;  
	}
	else 
		return std::vector<CSampleIDWeight>(0);
}

std::vector<CSampleIDWeight> CEquiEnergyModel::Initialize_MostDistant_WithinPercentile(const std::vector<CSampleIDWeight> &samples, int K, int stage_index, double percentile) const
{
	if ((int)samples.size() < K)
		return std::vector<CSampleIDWeight>(0); 	
	else if ((int)samples.size() == K)
		return samples; 
	
	CSampleIDWeight_Sorter comparator(parameter->lambda[stage_index]); 
	std::vector<CSampleIDWeight> sorted_samples(samples); 
	// descending on lambda*log_posterior
	std::sort(sorted_samples.begin(), sorted_samples.end(), comparator); 

	int n = (int)ceil((double)sorted_samples.size()*percentile); 
	if (n <= K)
		return std::vector<CSampleIDWeight>(sorted_samples.begin(), sorted_samples.begin()+K); 

	std::vector<CSampleIDWeight> starters(K); 
	// 1st starter, the best sample so far
	starters[0] = sorted_samples[0]; 
	// k-th starter, the most distant to starters(0..k-1) 
	double most_distance, distance; 
	TDenseVector e; 
	for (int ii=1; ii<K; ii++)
	{
		// examining all the n samples
		for(int jj=0; jj<n; jj++)
		{
			distance = 0.0; 
			// calculate the distance to starters[0, ..., ii-1]
			for (int kk=0; kk<ii; kk++)
			{
				e = samples[jj].data-starters[kk].data; 
				distance += InnerProduct(e,e); 
			}
			// use the sample that is fartherest to all existing centers
			if (jj == 0 || most_distance < distance)
			{
				most_distance = distance; 
				starters[ii] = sorted_samples[jj]; 
			}
		}
	}
	return starters; 
} 
