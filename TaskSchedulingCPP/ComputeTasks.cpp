#include <cstdlib>
#include <sstream>
#include <fstream>
#include "CEquiEnergyModel.h"
#include "CStorageHead.h"
#include "CEESParameter.h"
#include "CSampleIDWeight.h"
#include "dw_rand.h"

#include "slave_computing.h"

using namespace std; 

/*
    Normally, model.K(level+1) would have previously been set, but when restarting
    a job this is not the case.
*/
double ExecutingTuningTask(CEquiEnergyModel &model, int level, double Kplus)
{
  // setup
  model.Setup(level,Kplus);
	 
  // tune
  return model.Tune(true);
}

void ExecutingSimulationTask(CEquiEnergyModel &model, int node, double scale)
{
  // setup - in addition to fields set in ExecutingTuningTask()
  model.scale = scale;
  model.node = node;
  model.storage->SetupForOutput(model.energy_level,node);

  // simulate
  model.SimulateEE(true, model.parameter->Gn, model.parameter->N, false, true);

  // clean up
  model.storage->FlushOutput();

  // write diagnostic information
  model.WriteSimulationDiagnostic(model.energy_level, node);
}

void ExecutingInitialIndependentSimulationTask(CEquiEnergyModel &model, int node)
{
  // setup
  model.energy_level=model.parameter->number_energy_level;
  model.ReadInitialDistribution();
  model.storage->SetupForOutput(model.energy_level,node);

  // simulate
  int number_to_simulate = model.parameter->N * model.parameter->Gn;
  double log_kernel;
  TDenseVector x;
  TDraw y;
  for (int ii=0; ii < number_to_simulate; ii++)
    {
      log_kernel=model.InitialDraw(x);
      y.Set(x,model.target->LogPosterior(x),log_kernel,node,ii);
      model.storage->AddDraw(y); 
    }

  // clean up
  model.storage->FlushOutput();
}

double ExecutingInitialTuningTask(CEquiEnergyModel &model)
{
  // Setup
  model.energy_level=model.parameter->number_energy_level-1;
  model.scale=1.0;
  model.K(model.energy_level)=model.parameter->max_energy;
  model.storage->SetupForInput(model.energy_level+1,1,model.parameter->max_energy);
  model.storage->ComputeVarianceDecomposition(model.OrthonormalDirections,model.SqrtDiagonal);

  // tune
  return model.Tune(false);
}

void ExecutingInitialMetropolisSimulationTask(CEquiEnergyModel &model, double scale, int node)
{
  // setup - in addition to fields set in ExecutingInitialTuningTask()
  model.scale=scale;
  model.node=node;
  model.storage->SetupForOutput(model.energy_level,node);

  // burn-in and simulate
  int number_to_simulate = model.parameter->N * model.parameter->Gn;
  model.SimulateMH(false, number_to_simulate/5, false);
  model.SimulateMH(true, number_to_simulate, true);

  // clean up
  model.storage->FlushOutput();

  // write diagnostics
  // model.WriteSimulationDiagnostic(node);
}
