#include <cstdlib>
#include <sstream>
#include <fstream>
#include "CEquiEnergyModel.h"
#include "CStorageHead.h"
#include "CEESParameter.h"
#include "CSampleIDWeight.h"

#include "slave_computing.h"

using namespace std; 

// void ExecutingTuningTask(CEquiEnergyModel &model, int node)
// {
//   // prepare storage for obtaining draws from model.energy_level+1
//   model.storage->ClearStatus(model.energy_level); 
//   model.storage->ClearStatus(model.energy_level+1); 
//   model.storage->restore(model.energy_level);
//   model.storage->RestoreForFetch(model.energy_level+1);	

//   // Setup
//   cout << "calling model.ReadInitializationFile(model.energy_level);\n";
//   model.ReadInitializationFile(model.energy_level);
//   model.scale=1.0;
//   cout << model.SqrtDiagonal.dim << endl;

//   // tune
//   cout << "calling Tune(true)\n";
//   if (!model.Tune(true))
//     {
//       cerr << "Unable to tune MH sampler - level " << model.energy_level << endl;
//       abort();
//     }
//   cout << "returned from Tune(true)\n";

//   // write scale
//   model.WriteScale(model.energy_level,node,model.scale);

//   // // prepare storage
//   // model.storage->ClearStatus(model.energy_level);
//   // model.storage->ClearStatus(model.energy_level+1); 
//   // model.storage->finalize(model.energy_level);
//   // model.storage->ClearDepositDrawHistory(model.energy_level);
//   // model.storage->ClearDepositDrawHistory(model.energy_level+1);
// }

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

// void ExecutingSimulationTask(CEquiEnergyModel &model, int nNode, int node)
// {
//   // restore partial storage (previously obtained at this node) for updating
//   model.storage->ClearStatus(model.energy_level); 
//   model.storage->ClearStatus(model.energy_level+1); 
//   model.storage->restore(model.energy_level); 
//   model.storage->RestoreForFetch(model.energy_level+1);

//   // setup
//   model.ReadInitializationFile(model.energy_level);
//   int number_to_simulate = ceil(model.parameter->simulation_length/(nNode - 1));

//   // clean up
//   if (!model.SimulateEE(true, string(), number_to_simulate, model.parameter->reinitialize_factor, false, true))
//     {
//       cerr << "Error simulating model - level " << model.energy_level << endl;
//       abort();
//     }

//   // write diagnostics
//   model.WriteSimulationDiagnostic(node);

//   // finalze and clear-up storage
//   model.storage->ClearStatus(model.energy_level); 
//   model.storage->finalize(model.energy_level); 
//   model.storage->ClearDepositDrawHistory(model.energy_level);
//   model.storage->ClearDepositDrawHistory(model.energy_level+1); 
// }

void ExecutingInitialSimulationTask(CEquiEnergyModel &model, int nNode, int node)
{
  // Setup
  model.ReadInitializationFile(model.energy_level);
  int number_to_simulate = ceil(model.parameter->simulation_length/(nNode - 1));
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
