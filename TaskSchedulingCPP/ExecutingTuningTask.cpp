#include <cstdlib>
#include <sstream>
#include <fstream>
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
	convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.energy_level << "." << group_index; 
	string start_point_file = model.parameter->storage_dir + convert.str(); 
	if (!model.InitializeFromFile(start_point_file))
		return false; 

	// block scheme file
	convert.str(string()); 
	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_SCHEME; 
	string block_scheme_file_name = model.parameter->storage_dir  + convert.str(); 
	
	// tuning 
        convert.str(string());
       	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.energy_level << "." << group_index;
	string block_file_name = model.parameter->storage_dir + convert.str();	
	bool if_eejump=true; 
	model.storage->ClearStatus(model.energy_level); 
	model.storage->ClearStatus(model.energy_level+1); 
	model.storage->restore(model.energy_level);
	model.storage->RestoreForFetch(model.energy_level+1);	

	if (model.energy_level == model.parameter->number_energy_level-1)
	{
		if (!model.metropolis->AdaptiveBeforeSimulation_OnePass(model.current_sample, period, max_period, block_file_name, if_eejump, block_scheme_file_name))
		{
                	cerr << "CMetropolis::AdaptiveBeforeSimulation() : Error in writing " << block_file_name << endl;
                       	abort();
                }
	}
	else 
	{
		// based on samples of the previous level (not group specific) to estimate covariance matrix
		convert.str(string());
                convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level+1;
                string sample_file_name = model.parameter->storage_dir + convert.str();
                if (!model.metropolis->AdaptiveAfterSimulation_OnePass(model.current_sample, period, max_period, sample_file_name, block_file_name, if_eejump, block_scheme_file_name))
                {
                	cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl;
               	abort();
                }
	}
	return true; 
}


bool ExecutingTuningTask_AfterSimulation(size_t period, size_t max_period, CEquiEnergyModel &model, int group_index)
{
	// start point
	stringstream convert; 
	convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.energy_level << "." << group_index; 
	string start_point_file = model.parameter->storage_dir + convert.str(); 
	if (!model.InitializeFromFile(start_point_file))
		return false; 
	
	// block scheme file name
	convert.str(string());
        convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_SCHEME; 
        string block_scheme_file_name = model.parameter->storage_dir + convert.str();

	// block_file for writing
	convert.str(string());
       	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_2ND << model.energy_level << "." << group_index; 
        string block_file_name = model.parameter->storage_dir + convert.str();
	// sample file: based on the samples of the current level (group specific)
	convert.str(string());
        convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level << "." << group_index;	
	string sample_file_name = model.parameter->storage_dir + convert.str();
	bool if_eejump=true;
	model.storage->ClearStatus(model.energy_level); 
	model.storage->ClearStatus(model.energy_level+1); 
	model.storage->restore(model.energy_level);
	model.storage->RestoreForFetch(model.energy_level+1);

	if (!model.metropolis->AdaptiveAfterSimulation_OnePass(model.current_sample, period, max_period, sample_file_name, block_file_name, if_eejump, block_scheme_file_name) )
	{
		cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl;
		abort();
	}
	return true; 
}
