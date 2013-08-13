#include <cstdlib>
#include <vector>
#include <cmath>
#include <mpi.h>
#include "CEquiEnergy_TState.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

extern "C"
{
	#include "dw_parse_cmd.h"
}

using namespace std;

void master_deploying(int argc, char **argv, CEquiEnergy_TState &model, CEESParameter &parameter, CStorageHead &storage, const CSampleIDWeight &mode)
{	
	// 0: master node
	// 1:nNode-1: slave node
	size_t nNode;   
	MPI_Comm_size(MPI_COMM_WORLD, (int*)(&nNode)); 

	size_t n_initial = (size_t)dw_ParseInteger_String(argc, argv, "nInitial", nNode);
	
	vector<vector<unsigned int> > nodeGroup(n_initial);
        unsigned int i=1;
	while (i<nNode)
	{
        	for (unsigned int j=0; j<nodeGroup.size() && i<nNode; j++)
       	 	{
			nodeGroup[j].push_back(i);
                        i++;
                }
        }

	// make sure storage_directory is valid
	if (!storage.makedir())
	{
		cerr << "Error in making directory for " << parameter.run_id << endl; 
		double *sMessage= new double [N_MESSAGE];
       		for (unsigned int i=0; i<nodeGroup.size(); i++)
			for (unsigned int j=0; j<nodeGroup[i].size(); j++)
               			MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][j], END_TAG, MPI_COMM_WORLD);
       		delete [] sMessage;
		abort(); 
	}
	
	// Hill climb for the highest+1 level
	// number_hill_climb: number of times of hill climbing
	
	double received_log_posterior; 
	size_t number_hill_climb = (size_t) dw_ParseInteger_String(argc, argv, "HillClimb", 0); 
        if (number_hill_climb )
		received_log_posterior = DispatchHillClimbTask(nodeGroup, parameter, storage, number_hill_climb);

	// tuning and simulation
	int if_tuning_done = dw_ParseInteger_String(argc, argv, "TuningDone", 0);
	if (if_tuning_done == 0)	// if tuning is not done, need to tune
	{
		received_log_posterior = DispatchTuneSimulation(nodeGroup, parameter, storage, parameter.simulation_length, n_initial);  
	}
	else if (parameter.simulation_length > 0)
	{
		for (int level=parameter.highest_level; level>=parameter.lowest_level; level--)
		{
			cout << "Simulation at ... " << level << " ... for " << parameter.simulation_length << endl; 
			received_log_posterior = DispatchSimulation(nodeGroup, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG); 
	}
		}
	cout << "Done simulation" << endl; 

	// tell all the slaves to exit by sending an empty messag with 0 simulation length 
	double *sMessage= new double [N_MESSAGE];  
	for (unsigned int i=0; i<nodeGroup.size(); i++)
		for (unsigned int j=0; j<nodeGroup[i].size(); j++)
			MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][j], END_TAG, MPI_COMM_WORLD);
	delete [] sMessage; 
}

