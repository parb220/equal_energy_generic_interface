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
	double *sPackage= new double [N_MESSAGE];
	double *rPackage= new double [N_MESSAGE];   
	MPI_Status status;
	// Hill climb and GMM simulation at the highest+1 level
        if (nHillClimb )
	{
		DispatchHillClimbTask(nNode, model, nHillClimb);
		DispatchGMMSimulationTask(nNode, model, model.parameter->simulation_length); 
	}
	else // rebbing highest+1 level
	{
		model.storage->ClearStatus(model.parameter->highest_level+1);
        	// model.storage->consolidate(model.parameter->highest_level+1);
        	// model.storage->binning_equal_size(model.parameter->highest_level+1, model.parameter->number_energy_level, true, model.current_sample.GetSize_Data());
        	model.storage->binning_equal_size(model.parameter->highest_level+1, 50, true, model.current_sample.GetSize_Data());
        	model.storage->finalize(model.parameter->highest_level+1);
        	model.storage->ClearDepositDrawHistory(model.parameter->highest_level+1);

        	sPackage[LEVEL_INDEX] = (double)model.parameter->highest_level+1;
        	sPackage[RESERVE_INDEX] = (double)(model.storage->GetNumber_Bin(model.parameter->highest_level+1));
        	for (int i=0; i<model.storage->GetNumber_Bin(model.parameter->highest_level+1); i++)
                	sPackage[RESERVE_INDEX + i + 1] = model.storage->GetEnergyLowerBound(model.parameter->highest_level+1, i);

        	for (int i=1; i<nNode; i++)
                	MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, BINNING_INFO, MPI_COMM_WORLD);

        	for (int i=1; i<nNode; i++)
                	MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, BINNING_INFO, MPI_COMM_WORLD, &status);

	}

	// Simulation
	if (model.parameter->simulation_length) 
		DispatchTuneSimulation(nNode, nInitial, model, mode, model.parameter->simulation_length); 

	// tell all the slaves to exit by sending an empty messag with 0 simulation length 
	for (int i=1; i<nNode; i++)
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
        delete [] sPackage;
        delete [] rPackage;
}

