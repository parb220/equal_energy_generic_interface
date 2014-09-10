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

using namespace std; 

void DispatchHillClimbTask(int nNode, CEquiEnergyModel &model, int number_hill_climb)
{
	int nFeasibleSolutionPerNode = ceil((double)number_hill_climb/(double)nNode);

	double *sPackage = new double [N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = nFeasibleSolutionPerNode; 
	sPackage[LEVEL_INDEX] = model.parameter->number_energy_level ; 

	for (int i=1; i<nNode; i++)
	{
		sPackage[GROUP_INDEX] = i; 
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
	glob_t glob_result; 
	convert.str(string()); 
	convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE; 
	string filename = model.parameter->storage_dir + convert.str(); 
	
	glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result); 
	if (glob_result.gl_pathc == 1)
		rename(glob_result.gl_pathv[0], filename.c_str()); 
	else 
	{
		model.ClearGaussianMixtureModelParameters(); 
		for (int i=0; i<(int)glob_result.gl_pathc; i++)
		{
			if( !model.AggregateGaussianMixtureModelParameters(string(glob_result.gl_pathv[i])) )
			{
				cerr << "Error occurred while reading Gaussian mixture model parameters from " << glob_result.gl_pathv[i] << endl; 
				abort(); 
			}
			remove(glob_result.gl_pathv[i]);
		}
		if (!model.WriteGaussianMixtureModelParameters(filename))
		{
			cerr << "Error occurred while writing Gaussian mixture model parameters to " << filename << endl;
               		abort();
		}
	}
	globfree(&glob_result);

	// Consolidate partial storage files
	// model.storage->consolidate(model.parameter->number_energy_level);
}
