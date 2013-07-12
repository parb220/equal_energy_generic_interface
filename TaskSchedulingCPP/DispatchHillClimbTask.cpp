#include "master_deploying.h"

using namespace std; 

void DispatchHillClimbTask(const vector<unsigned int> &nodePool, const CEESParameter &parameter, CStorageHead &storage, int number_hill_climb)
{
	size_t nFeasibleSolutionPerNode = ceil((double)number_hill_climb/(double)(nodePool.size()));
        // Hill climber for the highest-plus-one level 

	double *sPackage = new double [N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = nFeasibleSolutionPerNode; 
	sPackage[BURN_INDEX] = 0; // irrelevant
	sPackage[FREQ_INDEX] = 0; // irrelevant
	sPackage[LEVEL_INDEX] = parameter.number_energy_level ; 
	sPackage[H0_INDEX] = parameter.h0; // irrelevant

	for (unsigned int i=0; i<nodePool.size(); i++)
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodePool[i], HILL_CLIMB_TAG, MPI_COMM_WORLD);		
	delete [] sPackage; 

	MPI_Status status; 
	double *rPackage = new double [N_MESSAGE];
	for (unsigned int i=0; i<nodePool.size(); i++)
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);
	delete [] rPackage;
 
	// Consolidate partial storage files
	unsigned int start_bin = parameter.BinIndex_Start(parameter.number_energy_level); 
	unsigned int end_bin = parameter.BinIndex_End(parameter.number_energy_level); 
	storage.consolidate(start_bin, end_bin);
}
