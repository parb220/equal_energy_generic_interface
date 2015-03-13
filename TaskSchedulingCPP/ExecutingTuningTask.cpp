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
	cout << "In ExecutingTuningTask_BeforeSimulation() - 1\n"; 
	// tuning 
        convert.str(string()); 
       	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.energy_stage << "." << group_index;
	string block_file_name = model.parameter->storage_dir + convert.str();	
	bool if_eejump;  
	model.storage->ClearStatus(model.energy_stage); 	
	model.storage->restore(model.energy_stage); 
	if (model.energy_stage < model.parameter->number_energy_stage)
	{
		model.storage->ClearStatus(model.energy_stage+1); cout << "In ExecutingTuningTask_BeforeSimulation() - 2\n";
		model.storage->RestoreForFetch(model.energy_stage+1);	cout << "In ExecutingTuningTask_BeforeSimulation() - 3\n";
	}	
	cout << "In ExecutingTuningTask_BeforeSimulation() - 4\n";
	if (model.energy_stage == model.parameter->number_energy_stage)
	{
		if_eejump=false; 
		if (!model.metropolis->AdaptiveBeforeSimulation_OnePass(model.current_sample, period, max_period, block_file_name, if_eejump, block_scheme_file_name))
		{
                	cerr << "ExecutingTuningTask_BeforeSimulation() : Error in writing " << block_file_name << endl;
                       	abort();
                }
	}
	else 
	{
		if_eejump= true;
		convert.str(string());
        	convert << model.parameter->run_id << "/" << model.parameter->run_id << BMATRIX;
		string bmatrix_file_name = model.parameter->storage_dir + convert.str(); 
		std::vector<TDenseMatrix> BMatrix = ReadBMatrixFile(bmatrix_file_name); 
		if (BMatrix.empty())
		{
			cerr << "ExecutingTuningTask_BeforeSimulation() : Error in reading " << bmatrix_file_name; 
			abort(); 
		}

		// if (!model.metropolis->AdaptiveAfterSimulation_WeightedSampling_OnePass(model.current_sample, period, max_period, samples, weight, block_file_name, if_eejump, block_scheme_file_name))
		if (!model.metropolis->AdaptiveAfterSimulation_WeightedSampling_OnePass(model.current_sample, period, max_period, BMatrix, block_file_name, if_eejump))
                {
                	cerr << "ExecutingTuningTask_BeforeSimulation() : Error in writing " << block_file_name << endl;
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

