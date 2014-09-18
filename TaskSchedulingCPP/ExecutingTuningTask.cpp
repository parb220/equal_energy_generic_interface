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
  // prepare storage for obtaining draws from model.energy_level+1
  model.storage->ClearStatus(model.energy_level); 
  model.storage->ClearStatus(model.energy_level+1); 
  model.storage->restore(model.energy_level);
  model.storage->RestoreForFetch(model.energy_level+1);	

  // Setup
  if (!model.SetupLevel(model.energy_level) || !model.SetScale(1.0))
    {
      cerr << "Unable to setup model for tuning for energy level " << model.energy_level << endl;
      abort(); 	
    }	

  // tune
  if (!model.Tune())
    {
      cerr << "Unable to tune MH sampler - level " << model.energy_level << endl;
      abort();
    }

  // write scale
  stringstream convert;
  convert << model.parameter->run_id << "/" << model.parameter->run_id << ".Scale." << model.energy_level << "." << group_index;
  string filename = model.parameter->storage_dir + convert.str();
  ofstream output_file;
  output_file.open(filename.c_str(), ios::out);
  if (!output_file)
    {
      cerr << "Error opening " << filename << endl; 
      abort(); 	
    }
  else
    {
      output_file << model.scale(model.energy_level);
    }
  output_file.close();

  // prepare storage
  model.storage->ClearStatus(model.energy_level);
  model.storage->ClearStatus(model.energy_level+1); 
  model.storage->finalize(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level+1);

  return true; 
}
