#include <vector>
#include <cmath>
#include <mpi.h>
#include <sstream>
#include "dw_rand.h"
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"

using namespace std; 

vector<string> glob(const string &pattern); 

void DispatchHillClimbTask(int nNode, int nInitial, CEquiEnergyModel &model, int number_hill_climb)
{
	int nFeasibleSolutionPerNode = ceil((double)number_hill_climb/(double)(nNode-1));

	double *sPackage = new double [N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = nFeasibleSolutionPerNode; 
	sPackage[LEVEL_INDEX] = model.parameter->number_energy_stage ; 

	for (int i=1; i<nNode; i++)
	{
		sPackage[GROUP_INDEX] = dw_uniform_int(nInitial);  
		sPackage[GROUP_NUMBER_INDEX] = 1; 
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, HILL_CLIMB_TAG, MPI_COMM_WORLD);		
	}
	delete [] sPackage; 

	MPI_Status status; 
	double *rPackage = new double [N_MESSAGE];
	for (int i=1; i<nNode; i++)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);
	}
	delete [] rPackage;
 
	// Consolidate gm_mean_hesssian files
	stringstream convert; 
	convert.str(string()); 
	convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE << ".*"; 
	string filename_pattern =  model.parameter->storage_dir + convert.str(); 
	convert.str(string()); 
	convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE; 
	string filename = model.parameter->storage_dir + convert.str(); 
	
	vector<string> merge_file_name = glob(filename_pattern); 
	if (merge_file_name.size() == 1)
		rename(merge_file_name[0].c_str(), filename.c_str()); 
	else 
	{
		model.ClearGaussianMixtureModelParameters(); 
		for (int i=0; i<(int)merge_file_name.size(); i++)
		{
			if( !model.AggregateGaussianMixtureModelParameters(merge_file_name[i]) )
			{
				cerr << "Error occurred while reading Gaussian mixture model parameters from " << merge_file_name[i] << endl; 
				abort(); 
			}
			remove(merge_file_name[i].c_str());
		}
		model.KeepOptimalGaussianMixtureModelParameters(); 
		if (!model.WriteGaussianMixtureModelParameters(filename))
		{
			cerr << "Error occurred while writing Gaussian mixture model parameters to " << filename << endl;
               		abort();
		}
	}
}
