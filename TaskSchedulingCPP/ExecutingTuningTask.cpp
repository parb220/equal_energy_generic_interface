#include <cstdlib>
#include <sstream>
#include <fstream>
#include "dw_math.h"
#include "CEquiEnergyModel.h"
#include "CStorageHead.h"
#include "CEESParameter.h"
#include "CMetropolis.h"
#include "CSampleIDWeight.h"
#include "storage_parameter.h"

using namespace std; 

bool ExecutingTuningTask_BeforeSimulation(size_t period, size_t max_period, CEquiEnergyModel &model, int group_index)
{
	// start point
	stringstream convert; 
	convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.energy_stage << "." << group_index; 
	string start_point_file = model.parameter->storage_dir + convert.str(); 
	if (!model.InitializeFromFile(start_point_file))
		return false; 

	// block scheme file
	convert.str(string()); 
	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_SCHEME; 
	string block_scheme_file_name = model.parameter->storage_dir  + convert.str(); 
	
	// tuning 
        convert.str(string());
       	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.energy_stage << "." << group_index;
	string block_file_name = model.parameter->storage_dir + convert.str();	
	bool if_eejump; 
	model.storage->ClearStatus(model.energy_stage); 
	model.storage->restore(model.energy_stage);
	if (model.energy_stage < model.parameter->number_energy_stage)
	{
		model.storage->ClearStatus(model.energy_stage+1); 
		model.storage->RestoreForFetch(model.energy_stage+1);	
	}	

	if (model.energy_stage == model.parameter->number_energy_stage)
	{
		if_eejump=false; 
		if (!model.metropolis->AdaptiveBeforeSimulation_OnePass(model.current_sample, period, max_period, block_file_name, if_eejump, block_scheme_file_name))
		{
                	cerr << "CMetropolis::AdaptiveBeforeSimulation() : Error in writing " << block_file_name << endl;
                       	abort();
                }
	}
	else 
	{
		if_eejump= true;
		// based on samples of the previous stage (not group specific) to estimate covariance matrix
		vector<CSampleIDWeight> samples; 
		if (!model.storage->DrawAllSample(model.energy_stage+1, samples) || samples.empty() )
		{
			cerr << "ExecutingTuningTask_BeforeSimulation() : Error in loading samples of previous stage.\n"; 
			abort(); 
		}
		vector<double> log_weight = model.Reweight(samples, model.energy_stage, model.energy_stage+1); 
		double log_weight_sum = log_weight[0]; 
		for (int i=1; i<(int)log_weight.size(); i++)
			log_weight_sum = AddLogs(log_weight_sum, log_weight[i]); 
		vector<double> weight(log_weight.size(), 0.0); 
		for (int i=0; i<(int)log_weight.size(); i++)
			weight[i] = exp(log_weight[i] - log_weight_sum); 

		if (!model.metropolis->AdaptiveAfterSimulation_WeightedSampling_OnePass(model.current_sample, period, max_period, samples, weight, block_file_name, if_eejump, block_scheme_file_name))
                {
                	cerr << "CMetroplis::AdaptiveAfterSimulation_WeightedSampling_OnePass() : Error in writing " << block_file_name << endl;
               		abort();
                }
	}

	model.storage->ClearStatus(model.energy_stage);
        model.storage->finalize(model.energy_stage);
        model.storage->ClearDepositDrawHistory(model.energy_stage);
	if (model.energy_stage < model.parameter->number_energy_stage)
	{
		model.storage->ClearStatus(model.energy_stage+1); 
        	model.storage->ClearDepositDrawHistory(model.energy_stage+1);
	}
	return true; 
}

