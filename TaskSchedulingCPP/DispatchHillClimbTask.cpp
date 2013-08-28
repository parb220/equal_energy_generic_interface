#include <vector>
#include <cmath>
#include <mpi.h>
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"

using namespace std; 

double DispatchHillClimbTask(const vector<vector<unsigned int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, int number_hill_climb)
{
	size_t nNode=0; 
	for (unsigned int i=0; i<nodeGroup.size(); i++)
		nNode += nodeGroup[i].size(); 

	size_t nFeasibleSolutionPerNode = ceil((double)number_hill_climb/(double)nNode);
	double max_log_posterior = -1.0e300, received_log_posterior; 

	double *sPackage = new double [N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = nFeasibleSolutionPerNode; 
	sPackage[BURN_INDEX] = 0; // irrelevant
	sPackage[FREQ_INDEX] = 0; // irrelevant
	sPackage[LEVEL_INDEX] = parameter.number_energy_level ; 
	sPackage[H0_INDEX] = parameter.h0; // irrelevant

	for (unsigned int i=0; i<nodeGroup.size(); i++)
	{
		for (unsigned int j=0; j<nodeGroup[i].size(); j++)
		{
			sPackage[GROUP_INDEX] = i; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][j], HILL_CLIMB_TAG, MPI_COMM_WORLD);		
		}
	}
	delete [] sPackage; 

	MPI_Status status; 
	double *rPackage = new double [N_MESSAGE];
	for (unsigned int i=0; i<nodeGroup.size(); i++)
	{
		for (unsigned int j=0; j<nodeGroup[i].size(); j++)
		{
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);
			received_log_posterior = rPackage[H0_INDEX]; 
			max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior;
		}
	}
	delete [] rPackage;
 
	// Consolidate partial storage files
	storage.consolidate(parameter.number_energy_level);
	return max_log_posterior; 
}
