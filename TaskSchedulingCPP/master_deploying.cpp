#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <mpi.h>
#include "CEquiEnergyModel.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

using namespace std;

void master_deploying(int nNode, int nHillClimb, int nInitial, CEquiEnergyModel &model, const CSampleIDWeight &mode)
{	
	// Hill climb and GMM simulation at the highest+1 level
        if (nHillClimb )
	{
	        DispatchHillClimbTask(nNode, model, nHillClimb);
	  	DispatchGMMSimulationTask(nNode, model, model.parameter->simulation_length); 
	}

	// Simulation
	if (model.parameter->simulation_length) 
	        DispatchTuneSimulation(nNode, nInitial, model, mode, model.parameter->simulation_length, false);   // original code last argument is not present

	// tell all the slaves to exit by sending an empty messag with 0 simulation length 
	double *sMessage= new double [N_MESSAGE];  
	for (int i=1; i<nNode; i++)
		MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
	delete [] sMessage; 
}

