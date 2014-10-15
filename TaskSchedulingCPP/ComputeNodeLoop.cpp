#include <mpi.h>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include "dw_dense_matrix.hpp"
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "mpi_parameter.h"

#include "ComputeNodeLoop.h"

using namespace std;

void ComputeNodeLoop(CEquiEnergyModel &model) 
{
  int my_rank, nNode; 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nNode);  
  MPI_Status status; 
	
  double *rPackage = new double [N_MESSAGE], *sPackage = new double [N_MESSAGE];    

  while (1)
    {
      MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
      if (status.MPI_TAG == END_TAG)
	{
	  cout << "Received END message" << endl;
	  delete []rPackage; 
	  delete []sPackage; 
	  exit(0); 
	}		
      else if (status.MPI_TAG == SIMULATION_TAG_INITIAL_INDEPENDENT)
	{
	  cout << "Received SIMULATION_TAG_INITIAL_INDEPENDENT message" << endl;

	  int node = (int)(rPackage[GROUP_INDEX]);

	  ExecutingInitialIndependentSimulationTask(model, node);
	}		
      else if (status.MPI_TAG == TUNE_TAG_INITIAL)
	{
	  cout << "Received TUNE_TAG_INITIAL message" << endl;

	  sPackage[SCALE_INDEX] = ExecutingInitialTuningTask(model);
	}		
      else if (status.MPI_TAG == SIMULATION_TAG_INITIAL_METROPOLIS)
	{
	  cout << "Received SIMULATION_TAG_INITIAL_METROPOLIS message" << endl;

	  int node = (int)(rPackage[GROUP_INDEX]);
	  double scale = rPackage[SCALE_INDEX];

	  ExecutingInitialMetropolisSimulationTask(model, scale, node);
	}
      else if (status.MPI_TAG == TUNE_TAG)
	{
	  cout << "Received TUNE_TAG message" << endl;

	  int level = (int)(rPackage[LEVEL_INDEX]);
	  double Kplus = rPackage[KPLUS_INDEX];

	  sPackage[SCALE_INDEX] = ExecutingTuningTask(model, level, Kplus);
	}
      else if (status.MPI_TAG == SIMULATION_TAG) 
	{	
	  cout << "Received SIMULATION_TAG message" << endl;

	  int node = (int)(rPackage[GROUP_INDEX]); 
	  double scale = rPackage[SCALE_INDEX];

	  ExecutingSimulationTask(model, node, scale); 
	}
      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
    }
  delete[] rPackage;
  delete[] sPackage;
}

/*
    Normally, model.K(level+1) would have previously been set, but when restarting
    a job this is not the case.
*/
double ExecutingTuningTask(CEquiEnergyModel &model, int stage)
{
  // setup
  model.SetupComputeNode(stage);
	 
  // tune
  return model.Tune(true);
}

void ExecutingSimulationTask(CEquiEnergyModel &model, double scale)
{
  // setup
  model.scale = scale;

  // simulate
  model.SimulateEE(true, model.parameter->Gn, model.parameter->N, false, true);

  // clean up
  model.storage->FinalizeComputeNode();

  // write diagnostic information
  model.WriteSimulationDiagnostic(model.stage, model.node);
}

void ExecutingInitialIndependentSimulationTask(CEquiEnergyModel &model)
{
  // setup
  model.SetupComputeNode(0);

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
  //model.storage->FlushOutput();
}

double ExecutingInitialTuningTask(CEquiEnergyModel &model)
{
  // Setup
  model.SetupComputeNode(1);

  // tune
  return model.Tune(false);
}

void ExecutingInitialMetropolisSimulationTask(CEquiEnergyModel &model, double scale)
{
  // setup - in addition to fields set in ExecutingInitialTuningTask()
  model.scale=scale;
  model.node=node;
  //model.storage->SetupForOutput(1,node);

  // burn-in and simulate
  int number_to_simulate = model.parameter->N * model.parameter->Gn;
  model.SimulateMH(false, number_to_simulate/5, false);
  model.SimulateMH(true, number_to_simulate, true);

  // clean up
  //model.storage->FlushOutput();

  // write diagnostics
  model.WriteSimulationDiagnostic(node);
}


