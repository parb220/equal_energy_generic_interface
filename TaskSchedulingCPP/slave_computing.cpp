#include <mpi.h>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include "dw_dense_matrix.hpp"
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CSampleIDWeight.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"

#include "slave_computing.h"

using namespace std;

void slave_computing(CEquiEnergyModel &model) 
{
  int my_rank, nNode; 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nNode);  
  MPI_Status status; 
	
  double *rPackage = new double [N_MESSAGE], *sPackage = new double [N_MESSAGE];    

  while (1)
    {
      cout << "Beginning poll loop" << endl;
      MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
      if (status.MPI_TAG == END_TAG)
	{
	  cout << "Received END message" << endl;
	  delete []rPackage; 
	  delete []sPackage; 
	  exit(0); 
	}		
      else if (status.MPI_TAG == BINNING_INFO)
	{
	  cout << "Received BINNING message" << endl;
	  model.energy_level = (int)(rPackage[LEVEL_INDEX]); 
	  int number_ring = (int)(rPackage[RESERVE_INDEX]); 
	  model.storage-> ResizeBin(model.energy_level, number_ring); 
	  for (int i=0; i<number_ring; i++)
	    model.storage->SetEnergyLowerBound(model.energy_level, i, rPackage[RESERVE_INDEX+i+1]);
	  if (model.energy_level+1 <= model.parameter->number_energy_level)
	    model.storage->ClearBin(model.energy_level+1); 
	  sPackage[RETURN_INDEX_1] = (double)(false); 
	}		
      else if (status.MPI_TAG == TUNE_TAG_INITIAL)
	{
	  cout << "Received TUNE_TAG_INITIAL message" << endl;
	  model.energy_level = (int)(rPackage[LEVEL_INDEX]);
	  model.node = (int)(rPackage[GROUP_INDEX]);

	  ExecutingInitialTuningTask(model, model.node);
	}		
      else if (status.MPI_TAG == SIMULATION_TAG_INITIAL)
	{
	  cout << "Received SIMULATION_TAG_INITIAL message" << endl;
	  model.energy_level = (int)(rPackage[LEVEL_INDEX]);
	  model.node = (int)(rPackage[GROUP_INDEX]);

	  ExecutingInitialSimulationTask(model, nNode, model.node);
	}
      else if (status.MPI_TAG == TUNE_TAG)
	{
	  cout << "Received TUNE message" << endl;
	  model.energy_level = (int)(rPackage[LEVEL_INDEX]);
	  model.node = (int)(rPackage[GROUP_INDEX]); 

	  ExecutingTuningTask(model, model.node);
	}
      else if (status.MPI_TAG == SIMULATION_TAG) 
	{	
	  cout << "Received SIMULATION message" << endl;
	  model.energy_level = (int)(rPackage[LEVEL_INDEX]);
	  model.node = (int)(rPackage[GROUP_INDEX]); 

	  ExecutingSimulationTask(model, nNode, model.node); 
	}
      cout << "Ending poll loop" << endl;
      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
    }
}

