#include <cmath>
#include <algorithm>
#include <vector> 
#include <functional>
#include <ctime>
#include <fstream>
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"

extern "C" {
	#include "dw_math.h"
	// #include "dw_switchio.h"
	#include "dw_rand.h"
}

using namespace std;

void CEquiEnergyModel::SaveSampleToStorage(CStorageHead &storage, const CSampleIDWeight &sample)
{
        storage.DepositSample(energy_level, storage.BinIndex(energy_level, -sample.weight), sample);
}

void CEquiEnergyModel::Take_Sample_Just_Drawn_From_Storage(const CSampleIDWeight &x)
{
	current_sample = x; 
	current_sample.id = (int)(time(NULL)-timer_when_started);
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

int CEquiEnergyModel::EE_Draw(const CEESParameter &parameter, CStorageHead &storage, size_t MH_thin)
{
	CSampleIDWeight x_new; 
	int new_sample_code = NO_JUMP; 
	double bounded_log_posterior_new; 

	if (energy_level == parameter.number_energy_level-1 || dw_uniform_rnd() > parameter.pee)
	{
		if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, MH_thin))
		{
			current_sample = x_new; 
			current_sample.id = (int)(time(NULL)-timer_when_started);
			new_sample_code = METROPOLIS_JUMP; 
		}
	}
	else 
	{
		// draw x_new from bin of the higher level of the same energy; 
		if (storage.DrawSample(energy_level+1, storage.BinIndex(energy_level+1,-current_sample.weight), x_new) ) // if a sample is successfully draw from bin
		{
			// calculate log_ratio in the current and the higher levels
			double log_ratio = parameter.LogRatio_Level(-x_new.weight, -current_sample.weight, energy_level); 
			log_ratio += parameter.LogRatio_Level(-current_sample.weight, -x_new.weight, energy_level+1); 
			if (log(dw_uniform_rnd()) <= log_ratio)
			{
				// accept the new sample
				Take_Sample_Just_Drawn_From_Storage(x_new); 
				new_sample_code = EQUI_ENERGY_JUMP; 
			}	
		}
		else	// when bin is empty, MH
		{
			if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, MH_thin) )
			{
				current_sample = x_new; 
				current_sample.id = (int)(time(NULL)-timer_when_started); 
				new_sample_code = METROPOLIS_JUMP; 
			}
		}
	}
	
	return new_sample_code; 
}

// A sample is randomly taken from a pool (with size desired_pool_size) of samples where the pool is formed by samples with higher log-posteriors. Note that samples with higher log-posterior values are stored in smaller-indexed bins. So, if the desired pool size is 10 while the size of the first bin is 100, then only the first bin will be used. In contrast, if the desired pool size is 100 while the total number of samples in the first 3 bins is barely greater than 100, then the first 3 bins will be used. 
bool CEquiEnergyModel::Initialize(CStorageHead &storage, size_t desiredPoolSize, int level_index)
{
	size_t N=storage.Number_Bin(level_index); 
	int nSample_Total=0;
	vector<int>nSample_Bin(N,0);
	for (int bin_index=0;  bin_index<N; bin_index++)
	{
		nSample_Bin[bin_index] = storage.GetNumberRecrod(level_index, bin_index); 
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
		if(storage.DrawSample(level_index, bin_index, current_sample))
		{
			current_sample.id = (int)(time(NULL)-timer_when_started); 
			// Because all samples stored in storage have had their log-posterior calculated and stored 
			// together with the sample values, there is no need to recalculate log-posterior at this moment
                	return true; 
		}	
	}	
	return false; 	
}

bool CEquiEnergyModel::InitializeWithBestSample(CStorageHead &storage, int level_index)
{
        int bin_index=0; 
        while (bin_index<storage.Number_Bin(level_index) && !(storage.DrawMostWeightSample(level_index, bin_index, current_sample) ) )
		bin_index ++; 
        if (bin_index >= storage.Number_Bin(level_index))
                return false;
	current_sample.id = (int)(time(NULL)-timer_when_started); 
        return true;
}

bool CEquiEnergyModel::InitializeWith_Kth_BestSample(CStorageHead &storage, size_t K, int level_index)
{
	vector<CSampleIDWeight> sample;
	int bin=0;  
	while (bin < storage.Number_Bin(level_index) && sample.size() < K)
	{
		storage.Draw_K_MostWeightSample(K, level_index, bin, sample); 
		bin ++; 
	}
	// if (sample.size() < K)
	//	return false; 
	Take_Sample_Just_Drawn_From_Storage(sample.back()); 
        return true;
}

bool CEquiEnergyModel::Initialize_RandomlyPickFrom_K_BestSample(CStorageHead &storage, size_t K, int level_index)
{
	vector<CSampleIDWeight> sample; 
	int bin = 0; 
	while (bin < storage.Number_Bin(level_index) && sample.size() < K)
	{
		storage.Draw_K_MostWeightSample(K, level_index, bin, sample); 
		bin ++; 
	}
	//if (sample.size() < K)
	//	return false; 
	int index = dw_uniform_int(K < sample.size() ? K : sample.size()); 
	Take_Sample_Just_Drawn_From_Storage(sample[index]); 
	return true; 
}


double CEquiEnergyModel::BurnIn(size_t burn_in_length)
{
	CSampleIDWeight x_new; 
	int nJump =0; 
	double max_posterior = current_sample.weight, bounded_log_posterior_new; 
	for (int i=0; i<burn_in_length; i++)
	{
		if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, 1) )
		{
			current_sample = x_new;
			current_sample.id = (int)(time(NULL)-timer_when_started);  
			max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior; 
			nJump ++; 
		}
	}
	cout << "MH Jump " << nJump << " out of " << burn_in_length << " in burning-in.\n"; 
	return max_posterior;  
}

double CEquiEnergyModel::Simulation_Within(const CEESParameter &parameter, CStorageHead &storage, bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new; 
	int nJump =0; 
	double max_posterior = current_sample.weight, bounded_log_posterior_new; 
	bool if_write_file = false; 
	ofstream output_file; 
	if (!sample_file_name.empty() )
	{
		output_file.open(sample_file_name.c_str(), ios::binary | ios::out); 
		if (output_file)
			if_write_file = true; 
	}

	for (int i=0; i<parameter.simulation_length; i++)
	{
		for (int j=0; j<parameter.thin; j++)
		{
			if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, parameter.THIN/parameter.thin) )
                	{
                        	current_sample = x_new;
                        	current_sample.id = (int)(time(NULL)-timer_when_started);
                        	max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior;
                        	nJump ++;
                	}
		}
		
		if (if_storage)
			SaveSampleToStorage(storage, current_sample);
		if (if_write_file)
			write(output_file, &current_sample); 
	}
	if (if_write_file)
		output_file.close(); 

	cout << "MH Jump " << nJump << " out of " << parameter.simulation_length*parameter.thin << " in simulation.\n"; 
	return max_posterior;  
}

double CEquiEnergyModel::Simulation_Cross(const CEESParameter &parameter, CStorageHead &storage, bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new;

	int nEEJump=0, nMHJump=0; 
	bool if_write_file = false;
        ofstream output_file;
        if (!sample_file_name.empty() )
        {
                output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
                if (output_file)
                        if_write_file = true;
        }

	double max_posterior = current_sample.weight; 
	for (int i=0; i<parameter.simulation_length; i++)
	{
		for (int j=0; j<parameter.thin; j++)
		{
			int jump_code = EE_Draw(parameter, storage, parameter.THIN/parameter.thin); 
			if (jump_code == EQUI_ENERGY_JUMP)
				nEEJump++; 
			else if (jump_code == METROPOLIS_JUMP)
				nMHJump++; 
			if (jump_code == EQUI_ENERGY_JUMP || jump_code == METROPOLIS_JUMP)
				max_posterior = max_posterior > current_sample.weight ? max_posterior : current_sample.weight; 
		}	
	
		if (if_storage)
			SaveSampleToStorage(storage, current_sample); 
		if (if_write_file)
			write(output_file, &current_sample); 
	}

	if (if_write_file)
		output_file.close(); 	

	cout << "EE Jump " << nEEJump << " out of " << parameter.simulation_length *parameter.thin<< " in simulation.\n"; 
	cout << "MH Jump " << nMHJump << " out of " << parameter.simulation_length *parameter.thin<< " in simulation.\n"; 
	return max_posterior; 
}

CEquiEnergyModel::CEquiEnergyModel() : 
if_bounded(true), energy_level(0), t_bound(1.0), current_sample(CSampleIDWeight()), metropolis(NULL), timer_when_started(time(NULL))
{}

CEquiEnergyModel::CEquiEnergyModel(bool _if_bounded, int eL, double _t, const CSampleIDWeight &_x, CMetropolis *_metropolis, time_t _time) :
if_bounded(_if_bounded), energy_level(eL), t_bound(_t), current_sample(_x), metropolis(_metropolis), timer_when_started(_time)
{
}

