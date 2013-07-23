#include <sstream>
#include "storage_parameter.h"
#include "master_deploying.h"

using namespace std;

double TopDownTuningSimulation(CEquiEnergyModel &model, const vector <unsigned int> &nodePool, const CEESParameter &parameter, CStorageHead &storage, const CSampleIDWeight &mode, size_t period, size_t max_period)
{
	string block_file_name, sample_file_name; 
	stringstream convert; 
	
	double max_log_posterior = -1.0e300, received_log_posterior; 

	for (unsigned int level = parameter.highest_level; level >= parameter.lowest_level; level --)
	{
		model.energy_level = level; 
		model.h_bound = parameter.h[level]; 
		model.t_bound = parameter.t[level]; 
		size_t estimation_length; 

		// storage of the previous (higher) level
		storage.RestoreForFetch(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) );
		// tuning is always based on the best sample so far.
		if (storage.empty(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) ) || !model.InitializeWithBestSample(storage, parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1)) )
                         model.current_sample = mode;
		storage.ClearDepositDrawHistory(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1)); 
	
		// file to save block information 	
		convert.str(string()); 
		convert << parameter.run_id << "/" << parameter.run_id << BLOCK_1ST << level; 
		block_file_name = parameter.storage_dir + convert.str(); 
		
		// STEP 1: 
		// if it is the highest temperature level, tune diag for the highest temp level
		// otherwise, tune  covariance based on samples at the higher temp level
		if (level == parameter.number_energy_level-1)
		{
			if (!model.metropolis->AdaptiveBeforeSimulation(model.current_sample, period, max_period, block_file_name))
			{
				cerr << "CMetropolis::AdaptiveBeforeSimulation() : Error in writing " << block_file_name << endl; 
				abort(); 
			}
		}
		else 
		{
			// file where samples are saved
			convert.str(string()); 
			convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << level+1; 
			sample_file_name = parameter.storage_dir + convert.str();  
			cout << "Tuning MH parameters (variance) based on the samples of the higher temperature level: " << level << endl; 
			if (!model.metropolis->AdaptiveAfterSimulation(model.current_sample, period, max_period, sample_file_name, block_file_name))
			{
				cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl; 
				abort(); 
			} 
		}

		// STEP 2: simulation to estimate covariance matrix of the current level 
		cout << "Simulating to estimate covariance matrix of the current level: " << level << endl; 
		estimation_length = 5000; 
		received_log_posterior = DispatchSimulation(nodePool, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_FIRST);	
		max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior; 

		
		// STEP 3: tuning covaraince based on the samples of the current level	
		convert.str(string()); 
		convert << parameter.run_id << "/" << parameter.run_id << BLOCK_2ND << level; 
		block_file_name = parameter.storage_dir + convert.str(); 
		convert.str(string());
                convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << level;
                sample_file_name = parameter.storage_dir + convert.str();
		cout << "Tuning MH parameters (variance) based on the samples of the current level : " << level << endl; 
		
		// start point for adaptive
		if (storage.empty(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) ) || !model.InitializeWithBestSample(storage, parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1)) )
                                model.current_sample = mode;
		if (!model.metropolis->AdaptiveAfterSimulation(model.current_sample, period, max_period, sample_file_name, block_file_name) )
		{
			cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl;
                        abort();
		} 

		// STEP 4: simulation
		cout << "Simulating for " << parameter.simulation_length << " steps at " << level << endl;
		recieved_log_posterior = DispatchSimulation(nodePool, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG); 
		max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior;

		// STEP 5: simulation to estimate covariance matrix of the lower temp level
		estimation_length = 5000; 
		cout << "Simulating to estimate covaraince matrix of the lower temp level : " << level << endl;
		received_log_posterior = DispatchSimulation(nodePool, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_SECOND);
		max_log_posterior = max_log_posterior < received_log_posterior ? max_log_posterior : received_log_posterior;
	}
	return max_log_posterior; 
}
