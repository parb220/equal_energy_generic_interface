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

  // setup
  // if (model.energy_level == model.parameter->number_energy_level-1)
  //   {
  //     model.energy_level++;
  //     model.ReadInitializationFile(model.energy_level);
  //     model.energy_level--;
  //   }
  model.ReadInitializationFile(model.energy_level);
  model.scale=1.0;

  // tune
  double scale = model.Tune();

  // write scale
  model.WriteScale(model.energy_level,group_index,scale);

  // prepare storage
  model.storage->ClearStatus(model.energy_level);
  model.storage->ClearStatus(model.energy_level+1); 
  model.storage->finalize(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level+1);

  return true; 
}
