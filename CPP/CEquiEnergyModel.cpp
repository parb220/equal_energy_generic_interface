#include <cmath>
#include <algorithm>
#include <vector> 
#include <functional>
#include <ctime>
#include "CEquiEnergyModel.h"

extern "C" {
	#include "dw_math.h"
	#include "dw_switchio.h"
	#include "dw_rand.h"
}

using namespace std;

void CEquiEnergyModel::SaveSampleToStorage(const CSampleIDWeight &sample, const CEESParameter &parameter, CStorageHead &storage)
{
	int binIndex = parameter.BinIndex(-sample.weight, energy_level);
        storage.DepositSample(binIndex, sample);
}

void CEquiEnergyModel::SaveSampleToFile(const CSampleIDWeight &sample, const CEESParameter &parameter, FILE *file)
{
	if (file == NULL)
                return;

	double log_prob = ((-sample.weight) > parameter.h[energy_level] ? (-sample.weight) : parameter.h[energy_level]) / parameter.t[energy_level]; 

        fprintf(file, "%le %le", log_prob, sample.weight);  
        for (int j=0; j<sample.data.dim; j++)
                fprintf(file, " %le", sample.data[j]);
        fprintf(file, "\n");
}

double CEquiEnergyModel::log_posterior_function(const double *x, int nx)
{
	double *old_x = new double[nx]; 
	// Save the old parameters stored in target_model to old_x
	ConvertThetaToFreeParameters(target_model, old_x); 
	ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) ); 
	
	// Post x to target_model
	ConvertFreeParametersToTheta(target_model, (double*)x); 
	ConvertFreeParametersToQ(target_model, (double*)x+NumberFreeParametersTheta(target_model) ); 
	double log_posterior = LogPosterior_StatesIntegratedOut(target_model);

	// Post old_x back to target_model 
	ConvertFreeParametersToTheta(target_model, old_x); 
	ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) ); 
	delete [] old_x; 

	return log_posterior;  
}

double CEquiEnergyModel::log_likelihood_function(const double *x, int nx)
{
	double *old_x = new double(nx); 
	// Save the old parameters stored in target_model to old_x
	ConvertThetaToFreeParameters(target_model, old_x); 
	ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) ); 

	// post x to target_model
	ConvertFreeParametersToTheta(target_model, (double*)x); 
	ConvertFreeParametersToQ(target_model, (double*)x+NumberFreeParametersTheta(target_model) ); 
	double log_likelihood = LogLikelihood_StatesIntegratedOut(target_model); 

	// post old_x back to target_model 
	ConvertFreeParametersToTheta(target_model, old_x); 
	ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) ); 
	delete [] old_x; 

	return log_likelihood; 
}

int CEquiEnergyModel::EE_Draw(const CEESParameter &parameter, CStorageHead &storage)
{
	CSampleIDWeight x_new; 
	int new_sample_code = NO_JUMP; 

	if (energy_level == parameter.number_energy_level-1 || dw_uniform_rnd() > parameter.pee)
	{
		if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data))
		{
			current_sample = x_new; 
			current_sample.id = (int)(time(NULL)-timer_when_started);
			new_sample_code = METROPOLIS_JUMP; 
		}
	}
	else 
	{
		// draw x_new from bin of the higher level of the same energy; 
		int bin_id = parameter.BinIndex(-current_sample.weight, energy_level+1);
		if (storage.DrawSample(bin_id, x_new)) // if a sample is successfully draw from bin
		{
			// calculate log_ratio in the current and the higher levels
			double log_ratio = parameter.LogRatio_Level(-x_new.weight, -current_sample.weight, energy_level); 
			log_ratio += parameter.LogRatio_Level(-current_sample.weight, -x_new.weight, energy_level+1); 
			if (log(dw_uniform_rnd()) <= log_ratio)
			{
				// accept the new sample
				current_sample = x_new; 
				current_sample.id = (int)(time(NULL)-timer_when_started);
				new_sample_code = EQUI_ENERGY_JUMP; 
			}	
			else 
			{
				if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data))
				{
					current_sample = x_new; 
					current_sample.id = (int)(time(NULL)-timer_when_started); 
					new_sample_code = METROPOLIS_JUMP; 
				}
			}
		}
		else 
		{
			if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data) )
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
bool CEquiEnergyModel::Initialize(CStorageHead &storage, unsigned int start_bin, unsigned int end_bin, size_t desiredPoolSize)
{
	size_t N=end_bin-start_bin+1; 
	unsigned int nSample_Total=0;
	vector<unsigned int>nSample_Bin(N,0);
	for (unsigned int bin=start_bin; bin<=end_bin; bin++)
	{
		nSample_Bin[bin-start_bin] = storage.GetNumberRecrod(bin); 
		nSample_Total += nSample_Bin[bin-start_bin]; 
	}

	unsigned int nLumSum = 0, bin= start_bin; 
	unsigned int random_index = dw_uniform_int(desiredPoolSize < nSample_Total ? desiredPoolSize : nSample_Total); 
	while (bin <= end_bin && !(random_index >= nLumSum && random_index < nLumSum + nSample_Bin[bin-start_bin]) )
	{	
		nLumSum += nSample_Bin[bin-start_bin]; 
		bin ++;
	}
	if (random_index >= nLumSum && random_index < nLumSum + nSample_Bin[bin-start_bin])
	{
		if(storage.DrawSample(bin, current_sample))
		{
			current_sample.id = (int)(time(NULL)-timer_when_started); 
			// Because all samples stored in storage have had their log-posterior calculated and stored 
			// together with the sample values, there is no need to recalculate log-posterior at this moment
			// current_sample.weight = log_posterior_function(current_sample.data); 
                	return true; 
		}	
	}	
	return false; 	
}

bool CEquiEnergyModel::InitializeWithBestSample(CStorageHead &storage, unsigned int start_bin, unsigned int end_bin)
{
        unsigned int bin=start_bin, number_sample;
        while (bin <= end_bin && (number_sample = storage.GetNumberRecrod(bin)) <= 0)
                bin ++;
        if (bin > end_bin)
                return false;
        if (storage.DrawLeastWeightSample(bin, current_sample))
	{
		current_sample.id = (int)(time(NULL)-timer_when_started); 
                return true;
	}
	else
        	return false;
}

double CEquiEnergyModel::BurnIn(size_t burn_in_length)
{
	CSampleIDWeight x_new; 
	unsigned int nJump =0; 
	double max_posterior = current_sample.weight; 
	for (unsigned int i=0; i<burn_in_length; i++)
	{
		if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data) )
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

double CEquiEnergyModel::Simulation_Within(const CEESParameter &parameter, CStorageHead &storage, bool if_storage)
{
	CSampleIDWeight x_new; 
	unsigned int nJump =0; 
	double max_posterior = current_sample.weight; 
	for (unsigned int i=0; i<parameter.simulation_length; i++)
	{
		for (unsigned int j=0; j<parameter.deposit_frequency; j++)
			if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data) )
                	{
                        	current_sample = x_new;
                        	current_sample.id = (int)(time(NULL)-timer_when_started);
                        	max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior;
                        	nJump ++;
                	}
		
		if (if_storage)
			SaveSampleToStorage(current_sample, parameter, storage); 
	}
	cout << "MH Jump " << nJump << " out of " << parameter.simulation_length << " in simulation.\n"; 
	return max_posterior;  
}

double CEquiEnergyModel::Simulation_Cross(const CEESParameter &parameter, CStorageHead &storage, bool if_storage)
{
	CSampleIDWeight x_new;

	unsigned int nEEJump=0, nMHJump=0; 

	double max_posterior = current_sample.weight; 
	for (unsigned int i=0; i<parameter.simulation_length; i++)
	{
		for (unsigned int j=0; j<parameter.deposit_frequency; j++)
		{
			int jump_code = EE_Draw(parameter, storage); 
			if (jump_code == EQUI_ENERGY_JUMP)
				nEEJump++; 
			else if (jump_code == METROPOLIS_JUMP)
				nMHJump++; 
			if (jump_code == EQUI_ENERGY_JUMP || jump_code == METROPOLIS_JUMP)
				max_posterior = max_posterior > current_sample.weight ? max_posterior : current_sample.weight; 
		}	
	
		if (if_storage)
			SaveSampleToStorage(current_sample, parameter, storage); 
	}
	
	cout << "EE Jump " << nEEJump << " out of " << parameter.simulation_length << " in simulation.\n"; 
	cout << "MH Jump " << nMHJump << " out of " << parameter.simulation_length << " in simulation.\n"; 
	return max_posterior; 
}

