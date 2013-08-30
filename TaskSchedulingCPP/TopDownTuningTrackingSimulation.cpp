#include <sstream>
#include <vector>
#include "CEquiEnergy_TState.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CSampleIDWeight.h"
#include "CMetropolis.h"
#include "storage_parameter.h"

using namespace std;

double TopDownTuningSimulation(CEquiEnergy_TState &model, const vector <vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, const CSampleIDWeight &mode, size_t period, size_t max_period)
{
	string block_file_name, sample_file_name; 
	stringstream convert; 
	
	double max_log_posterior = -1.0e300, received_log_posterior; 

	for (int level = parameter.highest_level; level >= parameter.lowest_level; level --)
	{
		model.energy_level = level; 
		model.t_bound = parameter.t[level]; 
		size_t estimation_length; 

		// storage of the previous (higher) level
		storage.RestoreForFetch(model.energy_level+1);
		// tuning is always based on the best sample so far.
		if (storage.empty(model.energy_level+1) || !model.InitializeWithBestSample(storage, model.energy_level+1) )
                         model.current_sample = mode;
		storage.ClearDepositDrawHistory(model.energy_level+1); 
	
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
		received_log_posterior = DispatchSimulation(nodeGroup, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_FIRST);	
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
		if (storage.empty(model.energy_level+1) || !model.InitializeWithBestSample(storage, model.energy_level+1) )
			model.current_sample = mode;
		if (!model.metropolis->AdaptiveAfterSimulation(model.current_sample, period, max_period, sample_file_name, block_file_name) )
		{
			cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl;
                        abort();
		} 

		// STEP 4: simulation
		cout << "Simulating for " << parameter.simulation_length << " steps at " << level << endl;
		recieved_log_posterior = DispatchSimulation(nodeGroup, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG); 
		max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior;

		// STEP 5: simulation to estimate covariance matrix of the lower temp level
		if (level > 0)
		{
			estimation_length = 5000; 
			cout << "Simulating to estimate covaraince matrix of the lower temp level : " << level << endl;
			received_log_posterior = DispatchSimulation(nodeGroup, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_SECOND);
			max_log_posterior = max_log_posterior < received_log_posterior ? max_log_posterior : received_log_posterior;
			storage.binning(level, parameter.number_energy_level, -log(0.5)/(1.0/parameter.t[level-1]-1.0/parameter.t[level]) ); 
			storage.finalize(level); 
			storage.ClearDepositDrawHistory(level);
		}
	}
	return max_log_posterior; 
}
