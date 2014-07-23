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
	vector<vector<int> > nodeGroup(nInitial);
        int i=1;
	while (i<nNode)
	{
        	for (int j=0; j<nodeGroup.size() && i<nNode; j++)
       	 	{
			nodeGroup[j].push_back(i);
                        i++;
                }
        }

	// Hill climb and GMM simulation at the highest+1 level
        if (nHillClimb )
	{
		DispatchHillClimbTask(nodeGroup, model, nHillClimb);
		DispatchGMMSimulationTask(nodeGroup, model, 1000*nHillClimb); 
	}

	// Simulation
	if (model.parameter->simulation_length) 
		DispatchTuneSimulation(nodeGroup, model, mode, model.parameter->simulation_length); 

	// tell all the slaves to exit by sending an empty messag with 0 simulation length 
	double *sMessage= new double [N_MESSAGE];  
	for (int i=0; i<nodeGroup.size(); i++)
		for (int j=0; j<nodeGroup[i].size(); j++)
			MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][j], END_TAG, MPI_COMM_WORLD);
	delete [] sMessage; 
}

