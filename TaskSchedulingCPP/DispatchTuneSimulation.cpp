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

using namespace std; 

void DispatchTuneSimulation(int nNode, CEquiEnergyModel &model, size_t simulation_length, bool save_space_flag)
{
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
	MPI_Status status; 

	if (model.parameter->highest_level >= model.parameter->number_energy_level)
	  {
	    cerr << "DispatchTuneSimulation(): model.parameter->highest_level must be less than model.parameter->number_energy_level" << endl;
	    abort();
	  }

	for (int level=model.parameter->highest_level; level>=model.parameter->lowest_level; level--)
	{
		sPackage[LEVEL_INDEX] = level;
		sPackage[PEE_INDEX] = model.parameter->pee; 

		model.storage->ClearStatus(level+1); 
		model.storage->RestoreForFetch(level+1); 

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		model.energy_level = level;
		model.SetupFromPreviousLevel(level);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Send out tuning jobs
		cout << "Dispatching tuning tasks to nodes - level " << level << endl;
		for (int i=1; i<nNode; i++)
		  {
		    sPackage[GROUP_INDEX] = i-1; 
		    MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG, MPI_COMM_WORLD);		
		  }
		for (int i=1; i<nNode; i++)
		  {
		    MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG, MPI_COMM_WORLD, &status);		
		  }

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Consolidate scales
		model.scale=model.ConsolidateScales(level);
		model.WriteInitializationFile();

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// simulation
		cout << "Simulation at " << level << endl; 
		sPackage[LENGTH_INDEX] = simulation_length;
	
		for (int i=1; i<nNode; i++)
		  {
		    sPackage[GROUP_INDEX] = i-1; 
		    MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG, MPI_COMM_WORLD);		
		  }	
		for (int i=1; i<nNode; i++)
		  {
		    MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG, MPI_COMM_WORLD, &status);
		  }

		model.WriteSimulationDiagnostic();

		// // write draws
		// if (level == model.parameter->number_energy_level-2)
		//   {
		//     fstream output_file;
		//     model.OpenFile(output_file,"Draws",level+1,true);
		//     vector<CSampleIDWeight> samples;  
		//     if (!model.storage->DrawAllSample(level+1, samples, true, model.current_sample.GetSize_Data()) || (samples.size() < model.parameter->simulation_length))
		//       {
		// 	cerr << "Error obtaining all samples from level " << level+1 << endl;
		// 	if (samples.size() < model.parameter->simulation_length) cerr << "Not enough samples - " << samples.size() << endl;
		// 	abort();
		//       }
		//     for (int ii=0; ii < samples.size(); ii++)
		//       output_file << samples[ii].weight << " " << samples[ii].data;
		//     output_file.close();
		//   }

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Consolidate partial storage files and bin
		model.storage->ClearStatus(level); 
		model.storage->consolidate(level); 
		model.storage->ClearStatus(level);
		model.storage->binning_equal_size(level, model.parameter->number_rings); 
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

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// to save space, remove level+1 samples
		if (save_space_flag && level+2 < model.parameter->highest_level-1 )
		{
	                stringstream convert;
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
