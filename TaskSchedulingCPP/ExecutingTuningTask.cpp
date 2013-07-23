#include "slave_computing.h"
#include "storage_parameter.h"

using namespace std; 

bool ExecutingTuningTask(size_t period, size_t max_period, CEquiEnergyModel &model, CStorageHead &storage, const CEESParameter &parameter, unsigned int my_rank, size_t initialPoolSize, const CSampleIDWeight &mode)
{
	storage.RestoreForFetch(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) );
	if (storage.empty(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) ) || !model.Initialize(storage, parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1), initialPoolSize) )
		model.current_sample = mode;
	storage.ClearDepositDrawHistory(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1));

	CSampleIDWeight tuning_start_point = model.current_sample; 

	string block_file_name, sample_file_name; 
	// block_file
	stringstream convert; 
	convert.str(string());
        convert << parameter.run_id << "/" << parameter.run_id << BLOCK_2ND << model.energy_level; 
        block_file_name = parameter.storage_dir + convert.str();
	// sample file
	convert.str(string());
        convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level+1;	
	sample_file_name = parameter.storage_dir + convert.str();

	size_t estimation_length = 5000; 
	size_t burn_in_length = 5000;  
	
	if (model.energy_level == parameter.number_energy_level-1)
	{
		if (!model.metropolis->FourPassAdaptive_StartWithoutSampleFile(tuning_start_point, period, max_period, estimation_length, burn_in_length, parameter.deposit_frequency, block_file_name))
			return false;  
	}
	else
	{ 
		if (!model.metropolis->FourPassAdaptive_StartWithSampleFile(tuning_start_point, period, max_period, estimation_length, burn_in_length, parameter.deposit_frequency, sample_file_name, block_file_name))
			return false;  
	}
	// tuning start point
	convert.str(string()); 
	convert << parameter.run_id << "/" << parameter.run_id << START_POINT << model.energy_level << "." << my_rank; 
	string start_point_file = parameter.storage_dir + convert.str();
	ofstream output_file;
	output_file.open(start_point_file.c_str(), ios::binary|ios::out); 
	if (!output_file)
		return false; 
	else 
		write(output_file, &tuning_start_point); 
	output_file.close(); 
	return true; 
}
