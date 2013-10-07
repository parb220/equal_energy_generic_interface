#include <cmath>
#include <vector>
#include <mpi.h>
#include <iostream>
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "storage_parameter.h"
#include "mpi_parameter.h"

void DispatchSimulation(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, int level, int tag);

using namespace std; 

void DispatchTuneSimulation(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length)
{
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
	MPI_Status status; 

	size_t estimation_length; 

	for (int level=parameter.highest_level; level>=parameter.lowest_level; level--)
	{
		sPackage[LEVEL_INDEX] = level; 
		sPackage[thin_INDEX] = parameter.thin; 
		sPackage[THIN_INDEX] = parameter.THIN; 
		sPackage[BURN_INDEX] = 0; 	// irrelevant
		sPackage[LENGTH_INDEX] = 0; 	// irrelevant

		// Tune before simulation 
		for (int i=0; i<nodeGroup.size(); i++)
		{
			sPackage[GROUP_INDEX] = i; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][0], TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD); 
		}
		for (int i=0; i<nodeGroup.size(); i++)
		{
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);
		}

		// Simulation to estimate group-specific covariance matrix
		estimation_length = 5000; 
		DispatchSimulation(nodeGroup, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_FIRST); 

		// Tune after simulation
		for (int i=0; i<nodeGroup.size(); i++)
		{
			sPackage[GROUP_INDEX] = i;
                        MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][0], TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD);
		}
		for (int i=0; i<nodeGroup.size(); i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD, &status);

		// simualtion
		cout << "Simulation at " << level << " for " << parameter.simulation_length << endl; 
		DispatchSimulation(nodeGroup, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG);
	
		// simualtion for the not-group-specific covariance of the lower temp level
		if (level > 0)
		{
			estimation_length = 5000;
			DispatchSimulation(nodeGroup, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_SECOND);

			// binning for the lower temperature-level's jump.  
			// storage.binning(level, parameter.number_energy_level, -log(0.5)/(1.0/parameter.t[level-1]-1.0/parameter.t[level])); 
			storage.binning_equal_size(level, parameter.number_energy_level); 
			storage.finalize(level); 
			storage.ClearDepositDrawHistory(level);
		}
	}

	delete []sPackage; 
	delete []rPackage; 
}
