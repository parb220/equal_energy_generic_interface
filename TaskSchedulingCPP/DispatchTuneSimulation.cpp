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

void DispatchSimulation(int nNode, int nInitial, CEquiEnergyModel &model, size_t simulation_length, int level, int message_tag); 

double EstimateLogMDD(CEquiEnergyModel &model, int level, int previous_level, double logMDD_previous);
double EstimateLogMDD(CEquiEnergyModel &model, int level, int proposal_type);

void DispatchTuneSimulation(int nNode, int nInitial, CEquiEnergyModel &model,const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag)
{
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
	MPI_Status status; 
	int iInitial;
	vector <int> availableNode;

	stringstream convert;
	string filename;
	ofstream output_file;
	ifstream input_file;

	size_t estimation_length; 

	for (int level=model.parameter->highest_level; level>=model.parameter->lowest_level; level--)
	{
		sPackage[LEVEL_INDEX] = level;
		sPackage[PEE_INDEX] = model.parameter->pee; 
		sPackage[thin_INDEX] = model.parameter->thin; 
		sPackage[THIN_INDEX] = model.parameter->THIN; 

		model.storage->ClearStatus(level+1); 
		model.storage->RestoreForFetch(level+1); 

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get independent directions from previous energy level draws
		if (!model.SetupLevel(level))
		  {
		    cerr << "Error computing independent directions from previous draws - energy level " << model.energy_level << endl;
		    abort();
		  }

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Send out tuning jobs
		cout << "Dispatching tuning tasks to nodes - level " << level << endl;
		availableNode.resize(nNode-1);
		for (int i=1; i<nNode; i++)
			availableNode[i-1] = i; 		
		iInitial = 0; 
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

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Consolidate scales
		double scale=0.0, scale_i;
		for (int i=0; i < nNode-1; i++)
		  {
		    convert.str(string());
		    convert << model.parameter->run_id << "/" << model.parameter->run_id << ".Scale." << level << "." << i;
		    filename = model.parameter->storage_dir + convert.str();
		    input_file.open(filename.c_str(), ios::in);
		    if (!input_file)
		      {
			cerr << "Error in opening " << filename << endl; 
			abort(); 	
		      }
		    else
		      {
			input_file >> scale_i;
			scale+=scale_i;
		      }
		    input_file.close();
		  }
		scale/=(double)(nNode-1);

		// write scale
      		convert.str(string());
		convert << model.parameter->run_id << "/" << model.parameter->run_id << ".Scale." << level;
		filename = model.parameter->storage_dir + convert.str();
		output_file.open(filename.c_str(), ios::out);
		if (!output_file)
		  {
		    cerr << "Error in opening " << filename << endl; 
		    abort(); 	
		  }
		else
		  output_file << scale;
		output_file.close();

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// simulation
		cout << "Simulation at " << level << " for " << model.parameter->simulation_length << endl; 
		sPackage[LENGTH_INDEX] = simulation_length;
		availableNode.resize(nNode-1);
		for (int i=1; i<nNode; i++)
			availableNode[i-1] = i; 
		iInitial = 0; 
		while (iInitial < nInitial)
		{
			if (availableNode.empty())
			{
				MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG, MPI_COMM_WORLD, &status);
				availableNode.push_back(status.MPI_SOURCE); 
			}
			// Assigns the last node in availableNode to iInitial, and remove the last node of availableNode
			sPackage[GROUP_INDEX] = iInitial; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, availableNode.back(), SIMULATION_TAG, MPI_COMM_WORLD); 
			availableNode.pop_back(); 
			iInitial ++; 
		}		
		for (int i=0; i<nNode-1-(int)availableNode.size(); i++)
		  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG, MPI_COMM_WORLD, &status);

		// Consolidate partial storage files and bin
		model.storage->ClearStatus(level); 
		model.storage->consolidate(level); 
		model.storage->ClearStatus(level);
		model.storage->binning_equal_size(level, 2*model.parameter->number_energy_level); 
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

		// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// // logMDD
	        // logMDD[level][level] = EstimateLogMDD(model, level, USE_TRUNCATED_POWER);
		// for (int j=level+1; j<(int)(logMDD[level].size()); j++)
		// 	logMDD[level][j] = EstimateLogMDD(model, level, level+1, logMDD[level+1][j]); 

		// cout << setprecision(20) << endl; 
		// cout << "logMDD at " << level << endl; 
		// for (int j=level; j<(int)(logMDD[level].size()); j++)
		// 	cout << logMDD[level][j] << "\t"; 
		// cout << endl; 

		// to save space, remove level+1 samples
		if (save_space_flag && level+2 < model.parameter->highest_level-1 )
		{
			model.storage->ClearSample(level+2);  
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << "*." << level+2 << ".*";;
                        string remove_file = model.parameter->storage_dir + convert.str();
			remove(remove_file.c_str()); 
		}

		if (model.K(model.energy_level) < 1.0000001) break;
	}

	delete []sPackage; 
	delete []rPackage; 
}
