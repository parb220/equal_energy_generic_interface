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

void DispatchHillClimbTask(const vector<vector<int> > &nodeGroup, CEquiEnergyModel &model, int number_hill_climb)
{
	size_t nNode=0; 
	for (int i=0; i<(int)nodeGroup.size(); i++)
		nNode += nodeGroup[i].size(); 

	size_t nFeasibleSolutionPerNode = ceil((double)number_hill_climb/(double)nNode);

	double *sPackage = new double [N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = nFeasibleSolutionPerNode; 
	sPackage[BURN_INDEX] = 0; // irrelevant
	sPackage[thin_INDEX] = 0; // irrelevant
	sPackage[THIN_INDEX] = 0; // irrelevant
	sPackage[LEVEL_INDEX] = model.parameter->number_energy_level ; 

	for (int i=0; i<(int)nodeGroup.size(); i++)
	{
		for (int j=0; j<(int)nodeGroup[i].size(); j++)
		{
			sPackage[GROUP_INDEX] = i; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][j], HILL_CLIMB_TAG, MPI_COMM_WORLD);		
		}
	}
	delete [] sPackage; 

	MPI_Status status; 
	double *rPackage = new double [N_MESSAGE];
	for (int i=0; i<(int)nodeGroup.size(); i++)
	{
		for (int j=0; j<(int)nodeGroup[i].size(); j++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);
	}
	delete [] rPackage;
 
	// Consolidate gm_mean_hesssian files
	stringstream convert; 
	convert.str(string()); 
	convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE << ".*"; 
	string filename_pattern =  model.parameter->storage_dir + convert.str(); 
	glob_t glob_result; 
	glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result); 
	std::vector<TDenseVector> gm_mean, gm_mean_per_file;
	std::vector<TDenseMatrix> gm_covariance_sqrt, gm_covariance_sqrt_per_file;
	for (int i=0; i<(int)glob_result.gl_pathc; i++)
	{
		if( !ReadGaussianMixtureModelParameters(string(glob_result.gl_pathv[i]), gm_mean_per_file, gm_covariance_sqrt_per_file) )
		{
			cerr << "Error occurred while reading Gaussian mixture model parameters from " << glob_result.gl_pathv[i] << endl; 
			abort(); 
		}
		gm_mean.insert(gm_mean.begin(), gm_mean_per_file.begin(), gm_mean_per_file.end()); 
		gm_covariance_sqrt.insert(gm_covariance_sqrt.begin(), gm_covariance_sqrt_per_file.begin(), gm_covariance_sqrt_per_file.end()); 
	}

	convert.str(string()); 
	convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE; 
	string filename = model.parameter->storage_dir + convert.str(); 
	if (!WriteGaussianMixtureModelParameters(filename, gm_mean, gm_covariance_sqrt))
	{
		cerr << "Error occurred while writing Gaussian mixture model parameters to " << filename << endl;
               	abort();
	}

	// Consolidate partial storage files
	// model.storage->consolidate(model.parameter->number_energy_level);
}
