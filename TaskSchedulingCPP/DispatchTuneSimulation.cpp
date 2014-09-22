#include <cmath>
#include <vector>
#include <mpi.h>
#include <glob.h>
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

void DispatchSimulation(int nNode, int nInitial, CEquiEnergyModel &model, size_t simulation_length, int level, int message_tag); 

double EstimateLogMDD(CEquiEnergyModel &model, int level, int previous_level, double logMDD_previous);
double EstimateLogMDD(CEquiEnergyModel &model, int level, int proposal_type);

void DispatchTuneSimulation(int nNode, int nInitial, CEquiEnergyModel &model,const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag)
{
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
	MPI_Status status; 

	size_t estimation_length; 

	// start point for tuning
	vector<CSampleIDWeight> start_points(nInitial); 
	vector<vector<double> > logMDD(model.parameter->number_energy_level, vector<double>(model.parameter->number_energy_level, 0.0)); 
	for (int level=model.parameter->highest_level; level>=model.parameter->lowest_level; level--)
	{
		sPackage[LEVEL_INDEX] = level; 
		sPackage[thin_INDEX] = model.parameter->thin; 
		sPackage[THIN_INDEX] = model.parameter->THIN; 
		sPackage[PEE_INDEX] = model.parameter->pee; 

		model.storage->ClearStatus(level+1); 
		model.storage->RestoreForFetch(level+1); 

		////////////////////////////////////////////////////////////////////////////////
		// Starting points
		if (model.storage->empty(level+1) || !model.Initialize_WeightedSampling(nInitial, level+1, start_points))
		// if (model.storage->empty(level+1) || !model.Initialize_MostDistant_WithinPercentile(nInitial, level+1, start_points, 0.30) ) 
		// !model.Initialize_KMeansClustering(nInitial, level+1, start_points) )
		{
			for (int i=0; i<nInitial; i++)
				start_points[i] = mode; 
		} 
		model.storage->ClearDepositDrawHistory(level+1); 

		stringstream convert; 
		string start_point_file; 
        	ofstream output_file;
		for (int i=0; i<nInitial; i++)
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
		}

		/////////////////////////////////////////////////////////////////////////////////
		// Send out tuning jobs
		vector <int> availableNode(nNode-1); 
		for (int i=1; i<nNode; i++)
			availableNode[i-1] = i; 
		
		int iInitial = 0; 
		while (iInitial < nInitial)
		{
			if (availableNode.empty())
			{
				MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);
				availableNode.push_back(status.MPI_SOURCE); 
			}
			// Assigns the last node in availableNode to iInitial, and remove the last node of availableNode
			sPackage[GROUP_INDEX] = iInitial; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, availableNode.back(), TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD); 
			availableNode.pop_back(); 
			iInitial ++; 
		}
		
		for (int i=0; i<nNode-1-(int)availableNode.size(); i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Simulation to estimate group-specific covariance matrix
		estimation_length = 5000; 
		DispatchSimulation(nNode, nInitial, model, estimation_length, level, TUNE_TAG_SIMULATION_FIRST); 

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Tune after simulation
		availableNode.clear(); 
		availableNode.resize(nNode-1); 
		for (int i=1; i<nNode; i++)
			availableNode[i-1] = i; 

		iInitial = 0; 
		while (iInitial < nInitial)
		{
			if(availableNode.empty())
			{
				MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD, &status);
				availableNode.push_back(status.MPI_SOURCE); 
			}
			// Assigns the last node in availableNode to iInitial, and remove the last node of availableNode

			sPackage[GROUP_INDEX] = iInitial;
                        MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, availableNode.back(), TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD);
			availableNode.pop_back(); 
			iInitial ++; 
		}
		for (int i=0; i<nNode-1-(int)availableNode.size(); i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_AFTER_SIMULATION, MPI_COMM_WORLD, &status);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// simualtion
		cout << "Simulation at " << level << " for " << model.parameter->simulation_length << endl; 
		DispatchSimulation(nNode, nInitial, model, model.parameter->simulation_length, level, SIMULATION_TAG);
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// simualtion for the not-group-specific covariance of the lower temp level
		estimation_length = 5000;
		DispatchSimulation(nNode, nInitial, model, estimation_length, level, TUNE_TAG_SIMULATION_SECOND);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// binning for the lower temperature-level's jump.  
		// storage.binning(level, parameter.number_energy_level, -log(0.5)/(1.0/parameter.t[level-1]-1.0/parameter.t[level])); 
		model.storage->ClearStatus(level); 
		model.storage->binning_equal_size(level, model.parameter->number_energy_level); 
		model.storage->finalize(level); 
		model.storage->ClearDepositDrawHistory(level);

		sPackage[LEVEL_INDEX] = (double)level; 
		sPackage[RESERVE_INDEX] = (double)(model.storage->GetNumber_Bin(level)); 
		for (int i=0; i<model.storage->GetNumber_Bin(level); i++)
			sPackage[RESERVE_INDEX + i + 1] = model.storage->GetEnergyLowerBound(level, i); 
		
		for (int i=1; i<nNode; i++)
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, BINNING_INFO, MPI_COMM_WORLD); 

		for (int i=1; i<nNode; i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, BINNING_INFO, MPI_COMM_WORLD, &status);
		
		// to save space, remove level+1 samples
		if (save_space_flag  && level > 0 && level+1 < model.parameter->number_energy_level-1 )
		{
			model.storage->ClearSample(level+1);  
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << "*." << level+1 << ".*";;
                        string remove_file = model.parameter->storage_dir + convert.str();
			glob_t glob_result1, glob_result2;
        		glob(remove_file.c_str(), GLOB_TILDE, NULL, &glob_result1);
                	for (int i=0; i<(int)glob_result1.gl_pathc; i++)
                        	remove(glob_result1.gl_pathv[i]);
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << "*.sample." << level+1;
			remove_file = model.parameter->storage_dir + convert.str();
                        glob(remove_file.c_str(), GLOB_TILDE, NULL, &glob_result2);
                        for (int i=0; i<(int)glob_result2.gl_pathc; i++)
                                remove(glob_result2.gl_pathv[i]);
                        globfree(&glob_result1);
                        globfree(&glob_result2);
			
		}
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// logMDD
	        logMDD[level][level] = EstimateLogMDD(model, level, USE_TRUNCATED_POWER);
		for (int j=level+1; j<(int)(logMDD[level].size()); j++)
			logMDD[level][j] = EstimateLogMDD(model, level, level+1, logMDD[level+1][j]); 

		cout << setprecision(20) << endl; 
		cout << "logMDD at " << level << endl; 
		for (int j=level; j<(int)(logMDD[level].size()); j++)
			cout << logMDD[level][j] << "\t"; 
		cout << endl; 

	}

	delete []sPackage; 
	delete []rPackage; 
}
