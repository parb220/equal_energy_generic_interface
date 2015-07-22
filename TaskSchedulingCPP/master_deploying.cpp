#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <mpi.h>
#include "CEquiEnergyModel.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

using namespace std;

void master_deploying(int nNode, int nInitial, CEquiEnergyModel &model, const CSampleIDWeight &mode, int nGroup_NSE)
{	
	double *sPackage= new double [N_MESSAGE];
	double *rPackage= new double [N_MESSAGE];   
	MPI_Status status;
	
	// Re-binning highest_stage+1 if continues from a previous run (such that highest_stage < number_energy_stage -1)

	if (model.parameter->highest_stage < model.parameter->number_energy_stage - 1)
	{
		model.storage->binning_equal_size(model.parameter->highest_stage+1, model.parameter->number_striation, model.parameter->lambda[model.parameter->highest_stage+1]); 

        	sPackage[LEVEL_INDEX] = (double)model.parameter->highest_stage+1;
        	sPackage[RESERVE_INDEX] = (double)(model.storage->GetNumber_Bin(model.parameter->highest_stage+1));
        	for (int i=0; i<model.storage->GetNumber_Bin(model.parameter->highest_stage+1); i++)
                	sPackage[RESERVE_INDEX + i + 1] = model.storage->GetEnergyLowerBound(model.parameter->highest_stage+1, i);

        	for (int i=1; i<nNode; i++)
                	MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, BINNING_INFO, MPI_COMM_WORLD);

        	for (int i=1; i<nNode; i++)
                	MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, BINNING_INFO, MPI_COMM_WORLD, &status);

	}

	// Tune & Simulation
	if (model.parameter->simulation_length) 
		DispatchTuneSimulation(nNode, nInitial, model, mode, model.parameter->simulation_length, nGroup_NSE); 

	// tell all the slaves to exit by sending an empty messag with 0 simulation length 
	for (int i=1; i<nNode; i++)
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
        delete [] sPackage;
        delete [] rPackage;
}

