#include <cmath>
#include <vector>
#include <mpi.h>
#include <iostream>
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "storage_parameter.h"
#include "mpi_parameter.h"

double DispatchSimulation(const vector<vector<unsigned int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, unsigned int level, int tag);

using namespace std; 

double DispatchTuneSimulation(const vector<vector<unsigned int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, size_t n_initial)
{
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
	MPI_Status status; 

	double max_log_posterior,  received_log_posterior; 
	size_t estimation_length; 

	for (unsigned int level=parameter.highest_level; level>=parameter.lowest_level; level--)
	{
		sPackage[LEVEL_INDEX] = level; 
		sPackage[H0_INDEX] = parameter.h0;
		sPackage[FREQ_INDEX] = parameter.deposit_frequency;
		sPackage[BURN_INDEX] = 0; 	// irrelevant
		sPackage[LENGTH_INDEX] = 0; 	// irrelevant

		// Tune before simulation 
		for (unsigned int i=0; i<nodeGroup.size(); i++)
		{
			sPackage[GROUP_INDEX] = i; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][0], TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD); 
		}
		for (unsigned int i=0; i<nodeGroup.size(); i++)
		{
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);
		}

		// Simulation to estimate group-specific covariance matrix
		estimation_length = 5000; 
		max_log_posterior = DispatchSimulation(nodeGroup, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_FIRST); 

		// Tune after simulation
		for (unsigned int i=0; i<nodeGroup.size(); i++)
		{
			sPackage[GROUP_INDEX] = i;
                        MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][0], TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD);
		}
		for (unsigned int i=0; i<nodeGroup.size(); i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD, &status);

		// simualtion
		cout << "Simulation at " << level << " for " << parameter.simulation_length << endl; 
		received_log_posterior = DispatchSimulation(nodeGroup, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG);
                max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior;
	
		// simualtion for the not-group-specific covariance of the lower temp level
		if (level > 0)
		{
			estimation_length = 5000;
			received_log_posterior = DispatchSimulation(nodeGroup, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_SECOND);
			max_log_posterior = max_log_posterior < received_log_posterior ? max_log_posterior : received_log_posterior;	

			// binning for the lower temperature-level's jump.  
			storage.binning(level, parameter.number_energy_level, -log(0.5)/(1.0/parameter.t[level-1]-1.0/parameter.t[level])); 
		}
	}

	delete []sPackage; 
	delete []rPackage; 
	return max_log_posterior; 
}
