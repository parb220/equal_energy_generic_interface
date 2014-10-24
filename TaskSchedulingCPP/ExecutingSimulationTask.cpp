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

void ExecutingSimulationTask(CEquiEnergyModel &model, int my_rank, int group_index, const CSampleIDWeight &mode, int message_tag)
{
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
	
	// start_point 
	stringstream convert;
        convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.energy_stage << "." << group_index;
        string start_point_file = model.parameter->storage_dir + convert.str();
        if (!model.InitializeFromFile(start_point_file) )
		model.current_sample = mode;

	// metropolis
        convert.str(string());
       	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.energy_stage; 
       	string block_file_name = model.parameter->storage_dir + convert.str();
       	if (!model.metropolis->ReadBlocks(block_file_name) )
       	{
       		cerr << "CMetropolis::ReadBlocks() : Error occurred while reading " << block_file_name << endl;
       		abort();
       	}

	// burn-in
	model.BurnIn(model.parameter->burn_in_length); 

	// whether to write dw output file
	if (message_tag == TUNE_TAG_SIMULATION_FIRST)
	{
		string sample_file_name; 
		convert.str(string()); 
		convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_stage << "." << group_index << "." << my_rank; 
		sample_file_name = model.parameter->storage_dir + convert.str(); 
		model.Simulation_Within(false, sample_file_name); 
	}
	else if (model.energy_stage == model.parameter->number_energy_stage)
		model.Simulation_Within(true, string()); 
	else 
		model.Simulation_Cross(true, string()); 

	// finalze and clear-up storage
	model.storage->ClearStatus(model.energy_stage); 
	model.storage->finalize(model.energy_stage); 
	if (model.energy_stage < model.parameter->number_energy_stage)
	{
		model.storage->ClearDepositDrawHistory(model.energy_stage);
		model.storage->ClearDepositDrawHistory(model.energy_stage+1); 
	}
}
