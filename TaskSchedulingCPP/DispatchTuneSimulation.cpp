#include <cmath>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "CSampleIDWeight.h"
#include "CEquiEnergyModel.h"
#include "storage_parameter.h"
#include "mpi_parameter.h"
#include "dw_matrix.h"
#include "mdd.hpp"
#include "mdd_function.h"

using namespace std; 

void DispatchSimulation(const vector<vector<int> > &nodeGroup, CEquiEnergyModel &model, size_t simulation_length, int level, int message_tag); 

double EstimateLogMDD(CEquiEnergyModel &model, int level, int previous_level, double logMDD_previous);
double EstimateLogMDD(CEquiEnergyModel &model, int level, int proposal_type);

void DispatchTuneSimulation(const vector<vector<int> > &nodeGroup, CEquiEnergyModel &model,const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag)
{
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
	MPI_Status status; 

	size_t estimation_length; 

	// start point for tuning
	vector<CSampleIDWeight> start_points(nodeGroup.size()); 
	vector<vector<double> > logMDD(model.parameter->number_energy_level, vector<double>(model.parameter->number_energy_level, 0.0)); 
	for (int level=model.parameter->highest_level; level>=model.parameter->lowest_level; level--)
	{
		sPackage[LEVEL_INDEX] = level; 
		sPackage[thin_INDEX] = model.parameter->thin; 
		sPackage[THIN_INDEX] = model.parameter->THIN; 
		sPackage[PEE_INDEX] = model.parameter->pee/10; 

		model.storage->ClearStatus(level+1); 
		model.storage->RestoreForFetch(level+1); 
		if (model.storage->empty(level+1) || !model.Initialize_MostDistant_WithinPercentile(nodeGroup.size(), level+1, start_points, 0.30) ) // !model.Initialize_KMeansClustering(nodeGroup.size(), level+1, start_points) )
		{
			for (int i=0; i<(int)(nodeGroup.size()); i++)
				start_points[i] = mode; 
		} 
		model.storage->ClearDepositDrawHistory(level+1); 

		stringstream convert; 
		string start_point_file; 
        	ofstream output_file;
		for (int i=0; i<nodeGroup.size(); i++)
		{
        		convert.str(string());
        		convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << level << "." << i;
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

			sPackage[GROUP_INDEX] = i; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][0], TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD); 
		}
		for (int i=0; i<nodeGroup.size(); i++)
		{
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);
		}

		// Simulation to estimate group-specific covariance matrix
		estimation_length = 5000; 
		DispatchSimulation(nodeGroup, model, estimation_length, level, TUNE_TAG_SIMULATION_FIRST); 

		// Tune after simulation
		sPackage[PEE_INDEX] = model.parameter->pee/10; 
		for (int i=0; i<nodeGroup.size(); i++)
		{
			sPackage[GROUP_INDEX] = i;
                        MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][0], TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD);
		}
		for (int i=0; i<nodeGroup.size(); i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD, &status);

		// simualtion
		cout << "Simulation at " << level << " for " << model.parameter->simulation_length << endl; 
		DispatchSimulation(nodeGroup, model, model.parameter->simulation_length, level, SIMULATION_TAG);
		
		// simualtion for the not-group-specific covariance of the lower temp level
		estimation_length = 5000;
		DispatchSimulation(nodeGroup, model, estimation_length, level, TUNE_TAG_SIMULATION_SECOND);

		// binning for the lower temperature-level's jump.  
		// storage.binning(level, parameter.number_energy_level, -log(0.5)/(1.0/parameter.t[level-1]-1.0/parameter.t[level])); 
		model.storage->ClearStatus(level); 
		model.storage->binning_equal_size(level, model.parameter->number_energy_level); 
		model.storage->finalize(level); 
		model.storage->ClearDepositDrawHistory(level);

		// logMDD
	        logMDD[level][level] = EstimateLogMDD(model, level, USE_TRUNCATED_POWER);
		for (int j=level+1; j<(int)(logMDD[level].size()); j++)
			logMDD[level][j] = EstimateLogMDD(model, level, level+1, logMDD[level+1][j]); 

		cout << setprecision(20) << endl; 
		cout << "logMDD at " << level << endl; 
		for (int j=level; j<(int)(logMDD[level].size()); j++)
			cout << logMDD[level][j] << "\t"; 
		cout << endl; 

		

		// to save space, remove level+1 samples
		if (save_space_flag && level+2 < model.parameter->highest_level-1 )
			model.storage->ClearSample(level+2);  
	}

	delete []sPackage; 
	delete []rPackage; 
}
