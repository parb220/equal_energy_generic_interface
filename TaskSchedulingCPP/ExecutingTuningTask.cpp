#include "slave_computing.h"
#include "storage_parameter.h"

using namespace std; 
bool ExecutingTuningTask_BeforeSimulation(size_t period, size_t max_period, CEquiEnergyModel &model, CStorageHead &storage, const CEESParameter &parameter, unsigned int group_index, size_t pool_size, const CSampleIDWeight &mode)
{
	// start point
	storage.RestoreForFetch(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) );
	if (storage.empty(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1) ) || !model.Initialize_RandomlyPickFrom_K_BestSample(pool_size, storage, parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1)) )
		model.current_sample = mode;
	storage.ClearDepositDrawHistory(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1));

	// save the start point
	stringstream convert;
	convert.str(string()); 
	convert << parameter.run_id << "/" << parameter.run_id << START_POINT << model.energy_level << "." << group_index; 
	string start_point_file = parameter.storage_dir + convert.str();
	ofstream output_file;
	output_file.open(start_point_file.c_str(), ios::binary|ios::out); 
	if (!output_file)
		return false; 
	else 
		write(output_file, &model.current_sample); 
	output_file.close(); 
	
	// tuning 
        convert.str(string());
       	convert << parameter.run_id << "/" << parameter.run_id << BLOCK_1ST << model.energy_level << "." << group_index;
	string block_file_name = parameter.storage_dir + convert.str();	
	
	if (model.energy_level == parameter.number_energy_level-1)
	{
		if (!model.metropolis->AdaptiveBeforeSimulation(model.current_sample, period, max_period, block_file_name))
		{
                	cerr << "CMetropolis::AdaptiveBeforeSimulation() : Error in writing " << block_file_name << endl;
                       	abort();
                }
	}
	else 
	{
		// based on samples of the previous level (not group specific) to estimate covariance matrix
		convert.str(string());
                convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level+1;
                string sample_file_name = parameter.storage_dir + convert.str();
                if (!model.metropolis->AdaptiveAfterSimulation(model.current_sample, period, max_period, sample_file_name, block_file_name))
                {
                	cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl;
               	abort();
                }
	}
	return true; 
}


bool ExecutingTuningTask_AfterSimulation(size_t period, size_t max_period, CEquiEnergyModel &model, const CEESParameter &parameter, unsigned int group_index, const CSampleIDWeight &mode)
{
	// start point
	stringstream convert; 
	convert << parameter.run_id << "/" << parameter.run_id << START_POINT << model.energy_level << "." << group_index; 
	string start_point_file = parameter.storage_dir + convert.str(); 
	if (!model.InitializeFromFile(start_point_file))
		model.current_sample = mode; 

	// block_file for reading : unnecessary
	/*convert.str(string()); 
	convert << parameter.run_id << "/" << parameter.run_id << BLOCK_1ST << model.energy_level << "." << group_index;
	string block_file_name = parameter.storage_dir + convert.str();
	if (!model.metropolis->ReadBlocks(block_file_name) )
	{
		cerr << "CMetropolis::ReadBlocks() : Error occurred while reading " << block_file_name << endl;
                abort();
	}*/

	// block_file for writing
	convert.str(string());
       	convert << parameter.run_id << "/" << parameter.run_id << BLOCK_2ND << model.energy_level << "." << group_index; 
        block_file_name = parameter.storage_dir + convert.str();
	// sample file: based on the samples of the current level (group specific)
	convert.str(string());
        convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level << "." << group_index;	
	string sample_file_name = parameter.storage_dir + convert.str();

	if (!model.metropolis->AdaptiveAfterSimulation(model.current_sample, period, max_period, sample_file_name, block_file_name) )
	{
		cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl;
		abort();
	}
	return true; 
}
