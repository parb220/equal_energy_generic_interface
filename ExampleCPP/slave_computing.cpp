#include <mpi.h>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include "dw_dense_matrix.hpp"
#include "CEquiEnergy_TState.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CSampleIDWeight.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"
#include "ReadWriteGaussianMixtureModelParameters.hpp"

extern "C"
{
        #include "dw_parse_cmd.h"
}
#include "slave_computing.h"

using namespace std;

void slave_computing(int argc, char **argv, CEquiEnergy_TState &model, const CSampleIDWeight &mode) 
{
	int my_rank, nCPU; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &nCPU);
	
	MPI_Status status; 
	
	double *rPackage = new double [N_MESSAGE]; 
	double *sPackage = new double [N_MESSAGE];    
	
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
		else if (status.MPI_TAG == HILL_CLIMB_TAG)
		{
			std::vector<TDenseVector> gm_mean; 
			std::vector<TDenseMatrix> gm_covariance_sqrt;  
			model.HillClimb_NPSOL(rPackage[LENGTH_INDEX], gm_mean, gm_covariance_sqrt);
			// Save gm_mean and gm_covariance_sqrt into files 
			stringstream convert;
			convert.str(string()); 
			convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE << "." << my_rank; 
			string gm_file = model.parameter->storage_dir + convert.str();  
			if (!WriteGaussianMixtureModelParameters(gm_file, gm_mean, gm_covariance_sqrt))
			{
				cerr << "Error occurred while writing Gaussian mixture model parameters to " << gm_file << endl;
                                abort();
			}
		}
		else if (status.MPI_TAG == GAUSSIAN_MIXTURE_MODEL_SIMULATION)
		{
			// model.energy_level = (int)(rPackage[LEVEL_INDEX]);
			// model.storage->restore(model.energy_level);
			// model.storage->finalize(model.energy_level);
        		// model.storage->ClearDepositDrawHistory(model.energy_level);
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
			
			size_t period = (size_t)dw_ParseInteger_String(argc, argv, "pr", 20);
                	size_t max_period = (size_t)dw_ParseInteger_String(argc, argv, "mpr", 16*period);
			size_t n_initial = (size_t)dw_ParseInteger_String(argc, argv, "nInitial", nCPU);

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
				if_within = true; 
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

			bool simulation_flag = ExecutingSimulationTask(if_within, if_write_file, if_storage, model, my_rank, group_index, 2*(size_t)nCPU, mode, status.MPI_TAG); 

			if (!simulation_flag)
			{
				cerr << "ExecutingSimulationTask: Error in simulation.\n"; 
				abort(); 
			}
		}
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
	}
}

