#include "slave_computing.h"
#include "storage_parameter.h"

extern "C"
{
	#include "dw_parse_cmd.h"
}

bool ExecutingSimulationTask(double &min_energy, bool if_within, bool if_write_sample_file, bool if_storage, CEquiEnergyModel &model, CStorageHead &storage, const CEESParameter &parameter, unsigned int my_rank, size_t initialPoolSize, const CSampleIDWeight &mode, int message_tag)
{
	// restore partial storage (previously obtained at this node) for updating
	storage.restore(parameter.BinIndex_Start(model.energy_level), parameter.BinIndex_End(model.energy_level));
	// Since the samples will be drawn from the higher level
	// the higher level needs to be restored for fetch (for partial record file)
	storage.RestoreForFetch(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) );
	// model::current_sample
	if (storage.empty(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) ) || !model.Initialize(storage, parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1), initialPoolSize) )
		model.current_sample = mode;

	// metropolis
	stringstream convert;
        convert.str(string());
        if (message_tag == TUNE_TAG_SIMULATION_FIRST)
        	convert << parameter.run_id << "/" << parameter.run_id << BLOCK_1ST << model.energy_level;
        else
        	convert << parameter.run_id << "/" << parameter.run_id << BLOCK_2ND << model.energy_level;
        
       	string block_file_name = parameter.storage_dir + convert.str();
       	if (!model.metropolis->ReadBlocks(block_file_name) )
       	{
       		cerr << "CMetropolis::ReadBlocks() : Error occurred while reading " << block_file_name << endl;
       		abort();
       	}

	double temp_min_energy; 
	// burn-in
	temp_min_energy = min_energy = -model.BurnIn(parameter.burn_in_length); 

	// whether to write dw output file
	string sample_file_name; 
	if (if_write_sample_file)
	{
		convert.str(string()); 
		convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level << "." << my_rank; 
		sample_file_name = parameter.storage_dir + convert.str(); 
	}
	else 
		sample_file_name = string(); 
	
	// simulation 
	if (if_within)
		temp_min_energy = -model.Simulation_Within(parameter, storage, if_storage, sample_file_name); 
	else
		temp_min_energy = -model.Simulation_Cross(parameter, storage, if_storage, sample_file_name); 

	// finalze and clear-up storage
	storage.finalize(parameter.BinIndex_Start(model.energy_level), parameter.BinIndex_End(model.energy_level)); 
	storage.ClearDepositDrawHistory(parameter.BinIndex_Start(model.energy_level), parameter.BinIndex_End(model.energy_level));
	storage.ClearDepositDrawHistory(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1)); 

	min_energy = min_energy <= temp_min_energy ? min_energy : temp_min_energy; 
	return true; 
}
