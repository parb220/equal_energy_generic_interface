#include <cmath>
#include <algorithm>
#include <vector> 
#include <functional>
#include <ctime>
#include "CEquiEnergyModel.h"
#include <fstream>

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

void CEquiEnergyModel::SaveSampleToFile(const CSampleIDWeight &sample, FILE *file)
{
	if (file == NULL)
                return;

	double log_likelihood = log_likelihood_function(sample); 

        fprintf(file, "%le %le", log_likelihood, sample.weight);  
        for (int j=0; j<sample.data.dim; j++)
                fprintf(file, " %le", sample.data[j]);
        fprintf(file, "\n");
}

bool CEquiEnergyModel::InitializeFromFile(const string &file_name)
{
	ifstream input_file; 
	input_file.open(file_name.c_str()); 
	if (!input_file)
		return false; 
	size_t dim; 
	double value; 

	input_file >> dim; 
	current_sample.data.Resize(dim); 
	for (unsigned int i=0; i<current_sample.data.dim; i++)
	{
		input_file >> value; 
		current_sample.data[i] = value; 
	}
	input_file.close(); 
	current_sample.id = (unsigned int)(time(NULL)-timer_when_started); 
	double bounded_log_posterior = log_posterior_function(current_sample); 
	return true;  
}

double CEquiEnergyModel::log_posterior_function(CSampleIDWeight &x)
{
	double *old_x = new double[x.data.dim]; 
	// Save the old parameters stored in target_model to old_x
	ConvertThetaToFreeParameters(target_model, old_x); 
	ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) ); 
	
	// Post x to target_model
	ConvertFreeParametersToTheta(target_model, x.data.vector); 
	ConvertFreeParametersToQ(target_model, x.data.vector+NumberFreeParametersTheta(target_model) ); 
	x.weight = LogPosterior_StatesIntegratedOut(target_model);
	double bounded_log_posterior; 
	if (if_bounded)
		bounded_log_posterior = -((-x.weight > h_bound) ? (-x.weight) : h_bound)/t_bound; 
	else 
		bounded_log_posterior = x.weight; 

	// Post old_x back to target_model 
	ConvertFreeParametersToTheta(target_model, old_x); 
	ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) ); 
	delete [] old_x; 

	return bounded_log_posterior;  
}

double CEquiEnergyModel::log_likelihood_function(const CSampleIDWeight &x)
{
	double *old_x = new double[x.data.dim]; 
	// Save the old parameters stored in target_model to old_x
	ConvertThetaToFreeParameters(target_model, old_x); 
	ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) ); 

	// post x to target_model
	ConvertFreeParametersToTheta(target_model, x.data.vector); 
	ConvertFreeParametersToQ(target_model, x.data.vector+NumberFreeParametersTheta(target_model) ); 
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
	double bounded_log_posterior_new; 

	if (energy_level == parameter.number_energy_level-1 || dw_uniform_rnd() > parameter.pee)
	{
		if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample))
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
				if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample))
				{
					current_sample = x_new; 
					current_sample.id = (int)(time(NULL)-timer_when_started); 
					new_sample_code = METROPOLIS_JUMP; 
				}
			}
		}
		else 
		{
			if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample) )
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
	double max_posterior = current_sample.weight, bounded_log_posterior_new; 
	for (unsigned int i=0; i<burn_in_length; i++)
	{
		if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample) )
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
	double max_posterior = current_sample.weight, bounded_log_posterior_new; 
	for (unsigned int i=0; i<parameter.simulation_length; i++)
	{
		for (unsigned int j=0; j<parameter.deposit_frequency; j++)
			if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample) )
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

CEquiEnergyModel::CEquiEnergyModel() : 
if_bounded(true), energy_level(0), h_bound(0.0), t_bound(1.0), current_sample(CSampleIDWeight()), target_model(NULL), metropolis(NULL), timer_when_started(time(NULL))
{}

CEquiEnergyModel::CEquiEnergyModel(bool _if_bounded, unsigned int eL, double _h, double _t, const CSampleIDWeight &_x, TStateModel *_model, CMetropolis *_metropolis, time_t _time) :
if_bounded(_if_bounded), energy_level(eL), h_bound(_h), t_bound(_t), current_sample(_x), target_model(_model), metropolis(_metropolis), timer_when_started(_time)
{
}

CEquiEnergyModel::~CEquiEnergyModel()
{
}
