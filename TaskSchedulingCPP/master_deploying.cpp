#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <mpi.h>
#include "CEquiEnergyModel.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

using namespace std;

void master_deploying(int nNode, CEquiEnergyModel &model)
{	
  if ((model.parameter->N > 0) && (model.parameter->G > 0))
    {
      // setup directories
      if (!model.storage->makedir())
	{
	  cerr << "Error in making directory for " << model.parameter->run_id << endl;
	}
      else
	{
	  // Initial distribution
	  if (model.parameter->highest_level >= model.parameter->number_energy_level)
	    {
	      DispatchInitialSimulation(nNode, model);
	      model.parameter->highest_level=model.parameter->number_energy_level-1;
	    }

	  // Simulation
	  DispatchTuneSimulation(nNode, model, model.parameter->N * model.parameter->G, false); 
	}
    }

  // tell all the slaves to exit by sending an empty messag with 0 simulation length 
  double *sMessage= new double [N_MESSAGE];  
  for (int i=1; i<nNode; i++)
    MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
  delete [] sMessage; 
}

