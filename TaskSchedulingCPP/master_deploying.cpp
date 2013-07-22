#include "master_deploying.h"
#include "mpi_parameter.h"

using namespace std;

void master_deploying(int argc, char **argv, CEquiEnergyModel &model, CEESParameter &parameter, CStorageHead &storage, const CSampleIDWeight &mode)
{	
	// 0: master node
	// 1:nNode-1: slave node
	size_t nNode;   
	MPI_Comm_size(MPI_COMM_WORLD, (int*)(&nNode)); 
	vector <unsigned int > nodePool(nNode-1); 
	for (unsigned int i=1; i<nNode; i++)
		nodePool[i-1] = i; 

	// make sure storage_directory is valid
	if (!storage.makedir())
	{
		cerr << "Error in making directory for " << parameter.run_id << endl; 
		double *sMessage= new double [N_MESSAGE];
       		for (int i=0; i<(int)(nodePool.size()); i++)
               		MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, nodePool[i], END_TAG, MPI_COMM_WORLD);
       		delete [] sMessage;
		abort(); 
	}
	
	// Hill climb for the highest+1 level
	// number_hill_climb: number of times of hill climbing
	
	size_t number_hill_climb = (size_t) dw_ParseInteger_String(argc, argv, "HillClimb", 0); 
	storage.RestoreForFetch(parameter.BinIndex_Start(parameter.number_energy_level), parameter.BinIndex_End(parameter.number_energy_level) ); 
        if (number_hill_climb && storage.empty(parameter.BinIndex_Start(parameter.number_energy_level), parameter.BinIndex_End(parameter.number_energy_level)) )
		DispatchHillClimbTask(nodePool, parameter, storage, number_hill_climb);
	storage.ClearDepositDrawHistory(parameter.BinIndex_Start(model.energy_level+1), parameter.BinIndex_End(model.energy_level+1));
	
	// tuning and simulation
	int if_tuning_done = dw_ParseInteger_String(argc, argv, "TuningDone", 0);
	if (if_tuning_done == 0)	// if tuning is not done, need to tune
	{
		// size_t random_block_tuning_number = (size_t)dw_ParseInteger_String(argc, argv, "RandomBlockTuning", 0);
		size_t period = (size_t)dw_ParseInteger_String(argc, argv, "pr", 20);
	        size_t max_period = (size_t)dw_ParseInteger_String(argc, argv, "mpr", 16*period);
		TopDownTuningSimulation(model, nodePool, parameter, storage, mode, period, max_period); 
	}
	else if (parameter.simulation_length > 0)
	{
		for (int level=parameter.highest_level; level>=parameter.lowest_level; level--)
		{
			cout << "Simulation at ... " << level << " ... for " << parameter.simulation_length << endl; 
			DispatchSimulation(nodePool, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG); 
	}
		}
	cout << "Done simulation" << endl; 

	// tell all the slaves to exit by sending an empty messag with 0 simulation length 
	double *sMessage= new double [N_MESSAGE];  
	for (int i=0; i<(int)nodePool.size(); i++)
		MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, nodePool[i], END_TAG, MPI_COMM_WORLD);
	delete [] sMessage; 
}

