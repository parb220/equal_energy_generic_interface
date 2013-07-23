#include "storage_parameter.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

using namespace std; 

double DispatchTuneSimulation(const vector<unsigned int> &nodePool, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length)
{
	double *sPackage = new double[N_MESSAGE]; 
	MPI_Status status; 
	double *rPackage = new double[N_MESSAGE]; 

	double max_log_posterior = -1.0e300, received_log_posterior; 
	for (unsigned int level=parameter.highest_level; level>=parameter.lowest_level; level--)
	{
		sPackage[LEVEL_INDEX] = level; 
		sPackage[H0_INDEX] = parameter.h0;
		sPackage[FREQ_INDEX] = parameter.deposit_frequency;
		sPackage[BURN_INDEX] = 0; 	// irrelevant
		sPackage[LENGTH_INDEX] = 0; 	// irrelevant

		cout << "Tuning at " << level << endl; 
		for (unsigned int i=0; i<nodePool.size(); i++)
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodePool[i], TUNE_TAG, MPI_COMM_WORLD); 
		for (unsigned int i=0; i<nodePool.size(); i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG, MPI_COMM_WORLD, &status);
		// simualtion
		cout << "Simulation at " << level << " for " << parameter.simulation_length << endl; 
		received_log_posterior = DispatchSimulation(nodePool, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG);
                max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior;
	
		// simualtion for the covariance of the lower temp level
		size_t estimation_length = 5000;
		cout << "Simulation for covariance matrix at" << level << " for " << estimation_length << endl; 
		received_log_posterior = DispatchSimulation(nodePool, parameter, storage, estimation_length, level, TUNE_TAG_SIMULATION_SECOND);
		max_log_posterior = max_log_posterior < received_log_posterior ? max_log_posterior : received_log_posterior;	
	}

	delete []sPackage; 
	delete []rPackage; 
	return max_log_posterior; 
}
