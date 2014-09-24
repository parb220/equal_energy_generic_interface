#include <cstdlib>
#include <sstream>
#include <fstream>
#include "CEquiEnergyModel.h"
#include "CStorageHead.h"
#include "CEESParameter.h"
#include "CSampleIDWeight.h"

#include "slave_computing.h"

using namespace std; 

void ExecutingTuningTask(CEquiEnergyModel &model, int node)
{
  // prepare storage for obtaining draws from model.energy_level+1
  model.storage->ClearStatus(model.energy_level); 
  model.storage->ClearStatus(model.energy_level+1); 
  model.storage->restore(model.energy_level);
  model.storage->RestoreForFetch(model.energy_level+1);	

  // setup
  model.ReadInitializationFile(model.energy_level);
  model.scale=1.0;

  // tune
  double scale = model.Tune();

  // write scale
  model.WriteScale(model.energy_level,node,scale);

  // prepare storage
  model.storage->ClearStatus(model.energy_level);
  model.storage->ClearStatus(model.energy_level+1); 
  model.storage->finalize(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level+1);
}

void ExecutingInitialTuningTask(CEquiEnergyModel &model, int node)
{
  // Setup
  model.ReadInitializationFile(model.energy_level);
  model.scale=1.0;

  // Draw initial starting point
  model.InitialDraw(model.current_sample.data);
  model.current_sample.weight=model.target->LogPosterior(model.current_sample.data);
  model.current_sample.calculated=true;

  // tune
  double scale = model.Tune();

  // write scale
  model.WriteScale(model.energy_level,node,scale);
}

void ExecutingSimulationTask(CEquiEnergyModel &model, int nNode, int node)
{
  // restore partial storage (previously obtained at this node) for updating
  model.storage->ClearStatus(model.energy_level); 
  model.storage->ClearStatus(model.energy_level+1); 
  model.storage->restore(model.energy_level); 
  model.storage->RestoreForFetch(model.energy_level+1);

  // simulate
  model.ReadInitializationFile(model.energy_level);
  model.SimulateEE(true, string(), model.parameter->G, model.parameter->N, false, true);
  model.WriteSimulationDiagnostic(node);

  // finalze and clear-up storage
  model.storage->ClearStatus(model.energy_level); 
  model.storage->finalize(model.energy_level); 
  model.storage->ClearDepositDrawHistory(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level+1); 
}

void ExecutingInitialSimulationTask(CEquiEnergyModel &model, int nNode, int node)
{
  // Setup
  model.ReadInitializationFile(model.energy_level);
  int number_to_simulate = model.parameter->N * model.parameter->G;
  string external_filename=model.MakeFilename("draws",model.energy_level,node);

  // Draw initial starting point
  model.InitialDraw(model.current_sample.data);
  model.current_sample.weight=model.target->LogPosterior(model.current_sample.data);
  model.current_sample.calculated=true;

  // Burn-in
  model.SimulateMH(false,string(),number_to_simulate/5, true);

  // simulate
  model.SimulateMH(false, external_filename, number_to_simulate, true);

  // write diagnostics
  // model.WriteSimulationDiagnostic(node);
}
