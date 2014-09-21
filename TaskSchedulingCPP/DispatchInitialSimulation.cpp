

#include <vector>
#include <cmath>
#include <mpi.h>
#include <glob.h>
#include <sstream>
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"

using namespace std; 

void DispatchInitialSimulation(int nNode, CEquiEnergyModel &model, int simulation_length)
{
  int nParameters=model.parameter->nParameters;
  int simulationLengthPerNode = ceil((double)simulation_length/(double)(nNode-1));
  model.energy_level=model.parameter->number_energy_level;

  MPI_Status status;
  double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];
  sPackage[LENGTH_INDEX] = (double)simulationLengthPerNode; 
  sPackage[LEVEL_INDEX] = (double)(model.energy_level);

  // setup file
  model.OrthonormalDirections=Identity(nParameters);
  model.SqrtDiagonal=Ones(nParameters);
 
  bool not_done=true;
  int iterations=0;
  do
    {
      // create initialization file
      model.CreateInitializationFile(model.energy_level,1.0,1.0,model.OrthonormalDirections,model.SqrtDiagonal);

      // send out tune command to compute nodes
      for (int i=1; i<nNode; i++)
	{
	  sPackage[GROUP_INDEX] = i-1; 
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG_INITIAL, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG, MPI_COMM_WORLD, &status);		
	}

      // consolidate scales and create initialization file
      model.CreateInitializationFile(model.energy_level,1.0,model.ConsolidateScales(model.energy_level),model.OrthonormalDirections,model.SqrtDiagonal);

      // send out simulate command to compute nodes
      for (int i=1; i<nNode; i++)
	{
	  sPackage[GROUP_INDEX] = i-1; 
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG_INITIAL, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG, MPI_COMM_WORLD, &status);
	}

      // read draws, set , compute ESS
      double ESS = model.AnalyzeInitialDraws(model.energy_level);
      cout << "ESS = " << ESS << endl;

      // exit if ESS is large enough or too many iterations have been performed
      if (ESS >= 1000)
	{
	  not_done=false;
	}
    }
  while(not_done);

  // simulate and save top level draws
  TDenseVector x;
  CSampleIDWeight y;
  for (int ii=0; ii < simulation_length; ii++)
    {
      model.InitialDraw(x);
      y.data=x;
      y.weight=model.target->LogPosterior(x);
      y.id=(int)(time(NULL)-model.timer_when_started);
      y.calculated=true;
      model.SaveSampleToStorage(y);  
    }

  // clean up
  model.storage->ClearStatus(model.energy_level);
  model.storage->binning_equal_size(model.energy_level, model.parameter->number_rings);
  model.storage->finalize(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level); 

  sPackage[LEVEL_INDEX] = (double)(model.energy_level);
  sPackage[RESERVE_INDEX] = (double)(model.storage->GetNumber_Bin(model.energy_level)); 
  for (int i=0; i < model.storage->GetNumber_Bin(model.energy_level); i++)
    sPackage[RESERVE_INDEX + i + 1] = model.storage->GetEnergyLowerBound(model.energy_level, i); 
		
  for (int i=1; i<nNode; i++)
    MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, BINNING_INFO, MPI_COMM_WORLD); 

  for (int i=1; i<nNode; i++)
    MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, BINNING_INFO, MPI_COMM_WORLD, &status);

  delete[]sPackage;
  delete[]rPackage;
}

