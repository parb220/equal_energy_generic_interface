#include <vector>
#include <cmath>
#include <mpi.h>
#include <sstream>
#include "dw_rand.h"
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"

using namespace std; 

void DispatchGMMSimulationTask(int nNode, int nInitial, CEquiEnergyModel &model, int simulation_length)
{
	int simulationLengthPerNode = ceil((double)simulation_length/(double)nNode);

	double *sPackage = new double [N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = simulationLengthPerNode; 
	sPackage[thin_INDEX] = model.parameter->thin;
        sPackage[THIN_INDEX] = model.parameter->THIN;
	sPackage[LEVEL_INDEX] = model.parameter->number_energy_stage ; 

	for (int i=1; i<nNode; i++)
	{
		sPackage[GROUP_INDEX] = dw_uniform_int(nInitial); 
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, GMM_SIMULATION_TAG, MPI_COMM_WORLD);		
	}

	MPI_Status status; 
	double *rPackage = new double [N_MESSAGE];
	for (int i=1; i<nNode; i++)
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, GMM_SIMULATION_TAG, MPI_COMM_WORLD, &status);
 
	// Consolidate partial storage files
	model.storage->ClearStatus(model.parameter->number_energy_stage);  
	model.storage->consolidate(model.parameter->number_energy_stage);
	model.storage->binning_equal_size(model.parameter->number_energy_stage,model.parameter->number_striation);
	model.storage->finalize(model.parameter->number_energy_stage); 
	model.storage->ClearDepositDrawHistory(model.parameter->number_energy_stage); 

	sPackage[LEVEL_INDEX] = (double)model.parameter->number_energy_stage;
        sPackage[RESERVE_INDEX] = (double)(model.storage->GetNumber_Bin(model.parameter->number_energy_stage));
        for (int i=0; i<model.storage->GetNumber_Bin(model.parameter->number_energy_stage); i++)
        	sPackage[RESERVE_INDEX + i + 1] = model.storage->GetEnergyLowerBound(model.parameter->number_energy_stage, i);

        for (int i=1; i<nNode; i++)
        	MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, BINNING_INFO, MPI_COMM_WORLD);

        for (int i=1; i<nNode; i++)
        	MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, BINNING_INFO, MPI_COMM_WORLD, &status);
	delete [] sPackage; 
	delete [] rPackage;
}

