#include <cmath>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "dw_rand.h"
#include "dw_matrix.h"
#include "CSampleIDWeight.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"
#include "storage_parameter.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

using namespace std; 

vector<string> glob(const string &pattern); 

std::vector<CSampleIDWeight> HighestPlus1Stage_Prior(int nNode, int nInitial, CEquiEnergyModel &model)
{
	return DispatchSimulation(nNode, nInitial, model, model.parameter->simulation_length,  model.parameter->number_energy_stage, SIMULATION_PRIOR_TAG);
}

std::vector<CSampleIDWeight> HighestPlus1Stage(int nNode, int nInitial, CEquiEnergyModel &model)
{
	// HillClimb to aquire an appropriate starting point
	DispatchHillClimbTask(nNode, nInitial, model, nInitial);
	
	stringstream convert; 
        convert.str(string());
        convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE;
        string gm_file = model.parameter->storage_dir + convert.str();
        if (!model.ReadGaussianMixtureModelParameters(gm_file) )
        {
        	cerr << "Error occurred while reading Gaussian mixture model parameters from " << gm_file << endl;
                abort();
        }
                
	vector<CSampleIDWeight> start_points(nInitial); 
        string start_point_file; 
	ofstream output_file;
        for (int i=0; i<nInitial; i++)
        {
        	start_points[i] = CSampleIDWeight(model.GetGMM_Mean(0)); 
                convert.str(string());
                convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.parameter->number_energy_stage << "." << i;
                start_point_file = model.parameter->storage_dir + convert.str();
                output_file.open(start_point_file.c_str(), ios::binary|ios::out);
                if (!output_file)
                {
                	cerr << "Error in writing to " << start_point_file << endl;
                        abort();
                }
                else 
                	write(output_file, &(start_points[i]));
                output_file.close();
       }

	// Adaptive and simulate o aquire covariance matrix	
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];
        MPI_Status status;
	sPackage[LEVEL_INDEX] = model.parameter->number_energy_stage;
        sPackage[thin_INDEX] = model.parameter->thin;
        sPackage[THIN_INDEX] = model.parameter->THIN;
        sPackage[PEE_INDEX] = model.parameter->pee/(model.parameter->THIN/model.parameter->thin);
       
	// Only need to run adaptive on each computing node once, because the results will be aggregated 
	for (int i=1; i<nNode; i++)
        {
		sPackage[GROUP_INDEX] = dw_uniform_int(nInitial);
		sPackage[GROUP_NUMBER_INDEX] = 1; 
                MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD);
        }

	for (int i=1; i<nNode; i++)
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);

	// Aggregate the scales obtained from computing nodes
	convert.str(string()); 
	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.parameter->number_energy_stage << ".*";  
	string block_file_pattern = model.parameter->storage_dir + convert.str();
	vector<string> block_file = glob(block_file_pattern); 

	convert.str(string());
        convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.parameter->number_energy_stage;
	string block_file_name = model.parameter->storage_dir + convert.str();
	if (!model.metropolis->AggregateBlocksAndRemoveFiles(block_file, block_file_name)) 
	{
		cerr << "Error in reading " << block_file_pattern << " or writing " << block_file_name << endl; 
		abort(); 
	}

	// Simulate to estimate covariance matrix
	// int estimation_length = 5000; 
	// DispatchSimulation(nNode, nInitial, model, estimation_length, model.parameter->number_energy_stage, TUNE_TAG_SIMULATION_FIRST);

	// Dispatch simulation
	return DispatchSimulation(nNode, nInitial, model, model.parameter->simulation_length,  model.parameter->number_energy_stage, SIMULATION_TAG);
	// DispatchGMMSimulationTask(nNode, nInitial, model, model.parameter->simulation_length); 
}
