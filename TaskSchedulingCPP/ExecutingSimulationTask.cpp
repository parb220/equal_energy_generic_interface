#include <cstdlib>
#include <sstream>
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"
#include "CStorageHead.h"
#include "CEESParameter.h"
#include "CSampleIDWeight.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"

using namespace std; 

double ExecutingSimulationTask(CEquiEnergyModel &model, int my_rank, int group_index, int nGroup, const CSampleIDWeight &mode, int message_tag)
{
	model.storage->InitializeBin(model.energy_stage, model.current_sample.GetSize_Data()); 
	// restore partial storage (previously obtained at this node) for updating
	model.storage->ClearStatus(model.energy_stage); 
	model.storage->restore(model.energy_stage); 
	// Since the samples will be drawn from the higher stage 
	// the higher stage needs to be restored for fetch (for partial record file)
	if (model.energy_stage < model.parameter->number_energy_stage)
	{
		model.storage->ClearStatus(model.energy_stage+1); 
		model.storage->RestoreForFetch(model.energy_stage+1);
	}

	double metropolis_accpt_rate = 0.0; 

	if (message_tag == SIMULATION_PRIOR_TAG)
		model.Simulation_Prior(true, string());	
	else 
	{
		// start_point 
		stringstream convert;
        	convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.energy_stage; // << "." << group_index;
        	string start_point_file = model.parameter->storage_dir + convert.str();
		std::vector<CSampleIDWeight> start_points; 
		if (!LoadSampleFromFile(start_point_file, start_points)) 
		{
			cerr << "ExecutingSimulationTask(): Error occurred reading start point files " << start_point_file << endl; 
			abort(); 
		}
        	// if (!model.InitializeFromFile(start_point_file, group_index) )
		
		// metropolis
        	convert.str(string());
       		convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.energy_stage; 
       		string block_file_name = model.parameter->storage_dir + convert.str();
       		if (!model.metropolis->ReadBlocks(block_file_name) )
       		{
       			cerr << "ExectuingSimulationTask: Error occurred while reading " << block_file_name << endl;
       			abort();
       		}

		int eGroup = group_index+nGroup <= (int)(start_points.size()) ? nGroup : (int)start_points.size()-group_index; 
		for (int iGroup = 0; iGroup < eGroup; iGroup++)
		{
			model.timer_when_started = group_index+iGroup; 
			model.current_sample = start_points[group_index+iGroup]; 

			// burn-in
			metropolis_accpt_rate = model.BurnIn(model.parameter->burn_in_length); 

			// whether to write dw output file
			if (message_tag == TUNE_TAG_SIMULATION_FIRST)
			{
				string sample_file_name; 
				convert.str(string()); 
				convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_stage << "." << group_index+iGroup << "." << my_rank; 
				sample_file_name = model.parameter->storage_dir + convert.str(); 
				model.Simulation_Within(false, sample_file_name); 
			}
			else if (message_tag == SIMULATION_TAG)
			{
				if (model.energy_stage == model.parameter->number_energy_stage)
					model.Simulation_Within(true, string());
				else
					model.Simulation_Cross(true, string()); 
			}
		}
	}
	// finalze and clear-up storage
	model.storage->ClearStatus(model.energy_stage); 
	model.storage->finalize(model.energy_stage); 
	model.storage->ClearDepositDrawHistory(model.energy_stage);
	if (model.energy_stage < model.parameter->number_energy_stage)
		model.storage->ClearDepositDrawHistory(model.energy_stage+1); 

	return metropolis_accpt_rate; 
}
