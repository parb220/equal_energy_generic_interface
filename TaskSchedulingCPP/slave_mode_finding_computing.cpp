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

using namespace std;

void slave_mode_finding_computing(int n_initial, CEquiEnergyModel &model, const CSampleIDWeight &mode, int optimizationN, int perturbationN, double perturbationS) 
{
	int my_rank, nNode; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &nNode);  
	MPI_Status status; 
	
	double *rPackage = new double [N_MESSAGE], *sPackage = new double [N_MESSAGE];    
	int group_index; 

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
			// start point
			stringstream convert;
			convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << "0." << group_index;
			string start_point_file = model.parameter->storage_dir + convert.str();
			if (!model.InitializeFromFile(start_point_file))
			{
				cerr << "Error occurred while reading start point file " << start_point_file << endl;
				abort(); 
			}
			//
			model.HillClimb_NPSOL(rPackage[LENGTH_INDEX], optimizationN, perturbationN, perturbationS, 1.0, model.current_sample.data); // model.parameter->t[model.parameter->number_energy_level]);
			// Save gm_mean and gm_covariance_sqrt into files 
			convert.str(string()); 
			convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE << "." << my_rank; 
			string gm_file = model.parameter->storage_dir + convert.str();  
			if (!model.WriteGaussianMixtureModelMeanAscii(gm_file))
			{
				cerr << "Error occurred while writing Gaussian mixture model parameters to " << gm_file << endl;
                                abort();
			}
			sPackage[RETURN_INDEX_1] = (double)(false); 
		}
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
	}
}

