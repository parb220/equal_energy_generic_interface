#include <mpi.h>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include "dw_dense_matrix.hpp"
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CSampleIDWeight.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"

#include "slave_computing.h"

using namespace std;

void slave_computing(int period, int max_period, int n_initial, CEquiEnergyModel &model, const CSampleIDWeight &mode, int optimizationN, int perturbationN, double perturbationS) 
{
	int my_rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
	MPI_Status status; 
	
	double *rPackage = new double [N_MESSAGE], *sPackage = new double [N_MESSAGE];    
	int group_index; 
	bool if_within, if_write_file, if_storage;   

	while (1)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
		if (status.MPI_TAG == END_TAG)
		{
			delete []rPackage; 
			delete []sPackage; 
			exit(0); 
		}
		else if (status.MPI_TAG == BINNING_INFO)
		{
			model.energy_level = (int)(rPackage[LEVEL_INDEX]); 
			int number_ring = (int)(rPackage[RESERVE_INDEX]); 
			model.storage-> ResizeBin(model.energy_level, number_ring); 
			for (int i=0; i<number_ring; i++)
				model.storage->SetEnergyLowerBound(model.energy_level, i, rPackage[RESERVE_INDEX+i+1]);
			if (model.energy_level+1 <= model.parameter->number_energy_level)
				model.storage->ClearBin(model.energy_level+1); 
			sPackage[RETURN_INDEX_1] = (double)(false); 
		}		
		else if (status.MPI_TAG == HILL_CLIMB_TAG)
		{
			model.HillClimb_NPSOL(rPackage[LENGTH_INDEX], optimizationN, perturbationN, perturbationS, model.parameter->t[model.parameter->number_energy_level]);
			// Save gm_mean and gm_covariance_sqrt into files 
			stringstream convert;
			convert.str(string()); 
			convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE << "." << my_rank; 
			string gm_file = model.parameter->storage_dir + convert.str();  
			if (!model.WriteGaussianMixtureModelParameters(gm_file))
			{
				cerr << "Error occurred while writing Gaussian mixture model parameters to " << gm_file << endl;
                                abort();
			}
			sPackage[RETURN_INDEX_1] = (double)(false); 
		}
		else if (status.MPI_TAG == GMM_SIMULATION_TAG)
		{
			model.energy_level = (int)(rPackage[LEVEL_INDEX]);
			group_index = (int)(rPackage[GROUP_INDEX]);
                        if (!GetCommunicationParameter(rPackage, N_MESSAGE, model.parameter))
                        {
                                cout << "GetCommunicationParameter() : Error occurred.\n";
                                abort();
                        }

			stringstream convert;
                        convert.str(string());
                        convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE; 
			string gm_file = model.parameter->storage_dir + convert.str(); 
			if (!model.ReadGaussianMixtureModelParameters(gm_file) )
			{
				cerr << "Error occurred while reading Gaussian mixture model parameters from " << gm_file << endl; 
				abort(); 
			}
			model.GMM_Simulation(rPackage[LENGTH_INDEX]); 
			model.storage->ClearStatus(model.energy_level); 
			model.storage->finalize(model.energy_level);
        		model.storage->ClearDepositDrawHistory(model.energy_level);
			sPackage[RETURN_INDEX_1] = (double)(false); 
		}
		else if (status.MPI_TAG == TUNE_TAG_BEFORE_SIMULATION || status.MPI_TAG == TUNE_TAG_AFTER_SIMULATION) 
		{
			model.energy_level = (int)(rPackage[LEVEL_INDEX]);
			group_index = (int)(rPackage[GROUP_INDEX]); 
			if (!GetCommunicationParameter(rPackage, N_MESSAGE, model.parameter))
			{
				cerr << "GetCommunicationParameter() : Error occurred.\n"; 
				abort(); 
			}
			model.t_bound = model.parameter->t[model.energy_level];
			
			if (status.MPI_TAG == TUNE_TAG_BEFORE_SIMULATION)
			{
				if (!ExecutingTuningTask_BeforeSimulation(period, max_period, model, group_index) )
				{
					cerr << "ExecutingTuningTask_BeforeSimulation() : Error occurred :: sample file reading or block_file writing or start_tune_point writing error.\n"; 
					abort(); 
				}
			}
			else if (status.MPI_TAG == TUNE_TAG_AFTER_SIMULATION)
			{
				if (!ExecutingTuningTask_AfterSimulation(period, max_period, model, group_index) )
                                {
                                        cerr << "ExecutingTuningTask_AfterSimulation() : Error occurred :: start_tune_point file reading or sample file reading or block_file writing error.\n";
                                        abort();
                                }
			}
		}
		else if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST || status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND || status.MPI_TAG == SIMULATION_TAG) 
		{	
			model.energy_level = (int)(rPackage[LEVEL_INDEX]);
			group_index = (int)(rPackage[GROUP_INDEX]); 
			if (!GetCommunicationParameter(rPackage, N_MESSAGE, model.parameter))
			{
				cout << "GetCommunicationParameter() : Error occurred.\n"; 
				abort(); 
			}
			model.t_bound = model.parameter->t[model.energy_level];

			if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST)
				if_within = false; 
			else 
				if_within = false; 
	
			if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST || status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND)
				if_write_file = true; 
			else 
				if_write_file = false; 

			if (status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND || status.MPI_TAG == SIMULATION_TAG)
				if_storage = true;
			else 
				if_storage = false; 	

			ExecutingSimulationTask(if_within, if_write_file, if_storage, model, my_rank, group_index, 100*n_initial, mode, status.MPI_TAG); 
		}
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
	}
}

