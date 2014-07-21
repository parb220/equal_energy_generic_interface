#include <vector>
#include <cmath>
#include <mpi.h>
#include <glob.h>
#include <sstream>
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"
#include "ReadWriteGaussianMixtureModelParameters.hpp"

using namespace std; 

void DispatchGMMSimulationTask(const vector<vector<int> > &nodeGroup, CEquiEnergyModel &model, int simulation_length)
{
	int nNode=0; 
	for (int i=0; i<(int)nodeGroup.size(); i++)
		nNode += (int)nodeGroup[i].size(); 

	int simulationLengthPerNode = ceil((double)simulation_length/(double)nNode);

	double *sPackage = new double [N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = simulationLengthPerNode; 
	sPackage[BURN_INDEX] = 0; // irrelevant
	sPackage[thin_INDEX] = 0; // irrelevant
	sPackage[THIN_INDEX] = 0; // irrelevant
	sPackage[LEVEL_INDEX] = model.parameter->number_energy_level ; 

	for (int i=0; i<(int)nodeGroup.size(); i++)
	{
		for (int j=0; j<(int)nodeGroup[i].size(); j++)
		{
			sPackage[GROUP_INDEX] = i; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][j], GMM_SIMULATION_TAG, MPI_COMM_WORLD);		
		}
	}
	delete [] sPackage; 

	MPI_Status status; 
	double *rPackage = new double [N_MESSAGE];
	for (int i=0; i<(int)nodeGroup.size(); i++)
	{
		for (int j=0; j<(int)nodeGroup[i].size(); j++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, GMM_SIMULATION_TAG, MPI_COMM_WORLD, &status);
	}
	delete [] rPackage;
 
	// Consolidate partial storage files
	model.storage->consolidate(model.parameter->number_energy_level);
}
