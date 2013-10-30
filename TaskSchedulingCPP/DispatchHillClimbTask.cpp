#include <vector>
#include <cmath>
#include <mpi.h>
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"

using namespace std; 

void DispatchHillClimbTask(const vector<vector<int> > &nodeGroup, CEquiEnergyModel &model, int number_hill_climb)
{
	size_t nNode=0; 
	for (int i=0; i<nodeGroup.size(); i++)
		nNode += nodeGroup[i].size(); 

	size_t nFeasibleSolutionPerNode = ceil((double)number_hill_climb/(double)nNode);

	double *sPackage = new double [N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = nFeasibleSolutionPerNode; 
	sPackage[BURN_INDEX] = 0; // irrelevant
	sPackage[thin_INDEX] = 0; // irrelevant
	sPackage[THIN_INDEX] = 0; // irrelevant
	sPackage[LEVEL_INDEX] = model.parameter->number_energy_level ; 

	for (int i=0; i<nodeGroup.size(); i++)
	{
		for (int j=0; j<nodeGroup[i].size(); j++)
		{
			sPackage[GROUP_INDEX] = i; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][j], HILL_CLIMB_TAG, MPI_COMM_WORLD);		
		}
	}
	delete [] sPackage; 

	MPI_Status status; 
	double *rPackage = new double [N_MESSAGE];
	for (int i=0; i<nodeGroup.size(); i++)
	{
		for (int j=0; j<nodeGroup[i].size(); j++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);
	}
	delete [] rPackage;
 
	// Consolidate partial storage files
	model.storage->consolidate(model.parameter->number_energy_level);
}
