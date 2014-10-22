#include <mpi.h>
#include <cstdlib>
#include <fstream>

#include "dw_dense_matrix.hpp"
#include "CEquiEnergyModel.h"

#include "mpi_parameter.h"
#include "ComputeNodes.h"

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

	  int stage = (int)(rPackage[STAGE_INDEX]);
	  ExecutingInitialIndependentSimulationTask(model, stage);
	}		
      else if (status.MPI_TAG == TUNE_TAG_INITIAL)
	{
	  cout << "Received TUNE_TAG_INITIAL message" << endl;

	  int stage = (int)(rPackage[STAGE_INDEX]);
	  sPackage[SCALE_INDEX] = ExecutingInitialTuningTask(model, stage);
	}		
      else if (status.MPI_TAG == SIMULATION_TAG_INITIAL_METROPOLIS)
	{
	  cout << "Received SIMULATION_TAG_INITIAL_METROPOLIS message" << endl;

	  int stage = (int)(rPackage[STAGE_INDEX]);
	  double scale = rPackage[SCALE_INDEX];
	  ExecutingInitialMetropolisSimulationTask(model, stage, scale);
	}
      else if (status.MPI_TAG == TUNE_TAG)
	{
	  cout << "Received TUNE_TAG message" << endl;

	  int stage = (int)(rPackage[STAGE_INDEX]);
	  sPackage[SCALE_INDEX] = ExecutingTuningTask(model, stage);
	}
      else if (status.MPI_TAG == SIMULATION_TAG) 
	{	
	  cout << "Received SIMULATION_TAG message" << endl;

	  int stage = (int)(rPackage[STAGE_INDEX]);
	  double scale = rPackage[SCALE_INDEX];
	  ExecutingSimulationTask(model, stage, scale); 
	}
      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
    }
  delete[] rPackage;
  delete[] sPackage;
}

double ExecutingInitialTuningTask(CEquiEnergyModel &model, int stage)
{
  // check
  if (stage != 1)
    {
      cerr << "ExecutingInitialTuningTask(): stage must be one" << endl;
      abort();
    }

  // setup
  model.SetupComputeNode(stage);

  // tune
  return model.Tune(false);
}

double ExecutingTuningTask(CEquiEnergyModel &model, int stage)
{
  // setup
  model.SetupComputeNode(stage);
	 
  // tune
  return model.Tune(true);
}

void ExecutingSimulationTask(CEquiEnergyModel &model, int stage, double scale)
{
  // check
  if (stage != model.stage)
    {
      cerr << "ExecutingTuningTask(): stage not equal to model.stage" << endl;
      abort();
    }

  // setup
  model.scale = scale;

  // simulate
  model.SimulateEE(true, model.parameter->Gn, model.parameter->N, false, true);

  // clean up
  model.storage->FinalizeComputeNode();

  // write diagnostic information
  model.WriteDiagnosticComputeNode();
}

void ExecutingInitialMetropolisSimulationTask(CEquiEnergyModel &model, int stage, double scale)
{
  // check
  if ((stage != model.stage) || (stage != 1))
    {
      cerr << "ExecutingTuningTask(): stage not equal to 1" << endl;
      abort();
    }

  // setup
  model.scale=scale;

  // burn-in and simulate
  int number_to_simulate = model.parameter->N * model.parameter->Gn;
  model.SimulateMH(false, number_to_simulate/5, false);
  model.SimulateMH(true, number_to_simulate, true);

  // clean up
  model.storage->FinalizeComputeNode();

  // write diagnostic information
  //model.WriteSimulationDiagnostic(model.stage, model.node);
}

void ExecutingInitialIndependentSimulationTask(CEquiEnergyModel &model, int stage)
{
  // check
  if (stage != 0)
    {
      cerr << "ExecutingInitialIndependentSimulationTask(): stage not equal to 0" << endl;
      abort();
    }

  // setup
  model.SetupComputeNode(stage);

  // simulate
  int number_to_simulate = model.parameter->N * model.parameter->Gn;
  double log_kernel;
  TDenseVector x;
  TDraw y;
  for (int ii=0; ii < number_to_simulate; ii++)
    {
      log_kernel=model.InitialDraw(x);
      y.Set(x,model.target->LogPrior(x),model.target->LogLikelihood(x),log_kernel,model.node,ii);
      model.storage->AddDraw(y); 
    }

  // clean up
  model.storage->FinalizeComputeNode();
}
