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

void ExecutingSimulationTask(bool if_within, bool if_write_sample_file, bool if_storage, CEquiEnergyModel &model, int my_rank, int group_index, size_t initialPoolSize, const CSampleIDWeight &mode, int message_tag)
{
	// restore partial storage (previously obtained at this node) for updating
	model.storage->ClearStatus(model.energy_level); 
	model.storage->ClearStatus(model.energy_level+1); 
	model.storage->restore(model.energy_level); 
	// Since the samples will be drawn from the higher level
	// the higher level needs to be restored for fetch (for partial record file)
	model.storage->RestoreForFetch(model.energy_level+1);
	
	// model::current_sample
	stringstream convert;
        convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.energy_level << "." << group_index;
        string start_point_file = model.parameter->storage_dir + convert.str();
        if (!model.InitializeFromFile(start_point_file) && (model.storage->empty(model.energy_level+1) || !model.Initialize_RandomlyPickFrom_K_BestSample(initialPoolSize, model.energy_level+1) ) )
		model.current_sample = mode;

	// metropolis
        convert.str(string());
        if (message_tag == TUNE_TAG_SIMULATION_FIRST)
        	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.energy_level << "." << group_index;
        else
        	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_2ND << model.energy_level << "." << group_index; 
        
       	string block_file_name = model.parameter->storage_dir + convert.str();
       	if (!model.metropolis->ReadBlocks(block_file_name) )
       	{
       		cerr << "CMetropolis::ReadBlocks() : Error occurred while reading " << block_file_name << endl;
       		abort();
       	}

	// burn-in
	model.BurnIn(model.parameter->burn_in_length); 

	// whether to write dw output file
	string sample_file_name; 
	if (if_write_sample_file)
	{
		convert.str(string()); 
		convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level << "." << group_index << "." << my_rank; 
		sample_file_name = model.parameter->storage_dir + convert.str(); 
	}
	else 
		sample_file_name = string(); 

	// simulation 
	if (if_within)
		model.Simulation_Within(if_storage, sample_file_name); 
	else
		model.Simulation_Cross(if_storage, sample_file_name); 

	// finalze and clear-up storage
	model.storage->ClearStatus(model.energy_level); 
	model.storage->finalize(model.energy_level); 
	model.storage->ClearDepositDrawHistory(model.energy_level);
	model.storage->ClearDepositDrawHistory(model.energy_level+1); 

}
