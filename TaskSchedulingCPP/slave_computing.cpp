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

