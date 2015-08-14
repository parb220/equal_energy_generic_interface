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

std::vector<int> ExecutingSimulationTask(CEquiEnergyModel &model, int my_rank, int group_index, int nGroup, const CSampleIDWeight &mode, int message_tag)
{
	model.storage->InitializeBin(model.energy_stage); 
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

	std::vector<int> nTotalJump(2,0), nOneTimeJump(2,0); // nJump[0] = nEEJump, nJump[1] = nMHJump
	if (message_tag == SIMULATION_PRIOR_TAG)
		nTotalJump = model.Simulation_Prior(true, string());	
	else 
	{
		// start_point 
        	string start_point_file = model.parameter->storage_dir + model.parameter->run_id + string("/") + model.parameter->run_id + START_POINT;
		std::vector<CSampleIDWeight> start_points = LoadSampleFromFile(start_point_file);  
		if (start_points.empty()) 
		{
			cerr << "ExecutingSimulationTask(): Error occurred reading start point files " << start_point_file << endl; 
			abort(); 
		}
		
		// metropolis
       		string block_file_name = model.parameter->storage_dir + model.parameter->run_id + string("/") + model.parameter->run_id + BLOCK;
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
			nOneTimeJump =  model.BurnIn(model.parameter->burn_in_length); 

			// whether to write dw output file
			if (message_tag == TUNE_TAG_SIMULATION_FIRST)
				nOneTimeJump = model.Simulation_Within(false, string()); 
			else if (message_tag == SCALE_MATRIX_FIT_TAG)
			{
				if (model.energy_stage == model.parameter->number_energy_stage)
					nOneTimeJump = model.Simulation_Within(false, string());
				else
					nOneTimeJump = model.Simulation_Cross(false, string()); 
			}
			else if (message_tag == SIMULATION_TAG)
			{
				if (model.energy_stage == model.parameter->number_energy_stage)
					nOneTimeJump = model.Simulation_Within(true, string());
				else
					nOneTimeJump = model.Simulation_Cross(true, string()); 
			}
			nTotalJump[0] += nOneTimeJump[0];
                        nTotalJump[1] += nOneTimeJump[1];
		}
	}
	// finalze and clear-up storage
	model.storage->ClearStatus(model.energy_stage); 
	model.storage->finalize(model.energy_stage); 
	model.storage->ClearDepositDrawHistory(model.energy_stage);
	if (model.energy_stage < model.parameter->number_energy_stage)
		model.storage->ClearDepositDrawHistory(model.energy_stage+1); 

	return nTotalJump; 
}
