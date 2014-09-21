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
  model.storage->RestoreForFetch(model.energy_level+1);

  // whether to write dw output file
  string sample_file_name; 
  if (if_write_sample_file)
    {
      stringstream convert;
      convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level << "." << group_index << "." << my_rank; 
      sample_file_name = model.parameter->storage_dir + convert.str(); 
    }
  else 
    sample_file_name = string(); 

  // simulate
  model.ReadInitializationFile(model.energy_level);
  // if (!model.SetupLevel(model.energy_level) || !model.SetScale(-1.0))
  //   {
  //     cerr << "Unable to setup model for simulation from energy level " << model.energy_level << endl;
  //     abort(); 	
  //   }  

  int number_to_simulate = 1 + model.parameter->simulation_length/model.parameter->n_compute_cores;

  model.Simulate(if_storage, sample_file_name, number_to_simulate, 2000, false, true);

  model.WriteSimulationDiagnostic(group_index);

  // finalze and clear-up storage
  model.storage->ClearStatus(model.energy_level); 
  model.storage->finalize(model.energy_level); 
  model.storage->ClearDepositDrawHistory(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level+1); 
}
