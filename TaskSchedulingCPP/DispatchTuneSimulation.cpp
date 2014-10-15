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

void DispatchTuneSimulation(int nNode, CEquiEnergyModel &model, bool save_space_flag)
{
  double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
  MPI_Status status; 

  // If model.parameter->highest_level is less than  model.parameter->number_energy_level-1, 
  // then K(level+1) will not be properly set. This is a bug if a warm restart is performed.

  double lambda=0.0;
  for (int stage=1; lambda < 1.0; stage++)
    {
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Send out tuning jobs
      cout << "Dispatching tuning tasks to nodes - stage " << stage << endl;
      double scale=0.0;

      for (int i=1; i<nNode; i++)
	{
	  sPackage[STAGE_INDEX]=(double)stage;
	  sPackage[LAMBDA_MIN_INDEX] = model.lambda(stage);
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG, MPI_COMM_WORLD);		
	}

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Setup master while compute nodes tune
      //model.Setup(level, model.K(level+1));
      //model.ESS(level) = model.storage->ComputeEffectiveSampleSize();

      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG, MPI_COMM_WORLD, &status);	
	  scale+=rPackage[SCALE_INDEX];
	}
      model.scale=scale/=(double)(nNode - 1);

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // simulation
      cout << "Dispatching simulation tasks to nodes - stage " << stage << endl; 
	
      for (int i=1; i<nNode; i++)
	{
	  sPackage[GROUP_INDEX] = (double)(i-1);
	  sPackage[SCALE_INDEX] = scale;
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG, MPI_COMM_WORLD);		
	}	
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG, MPI_COMM_WORLD, &status);
	}
      //      model.storage->Consolidate(level);
      //      model.WriteSimulationDiagnostic(level);

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // // to save space, remove level+1 samples
      // if (save_space_flag && level+2 < model.parameter->highest_level-1 )
      // {
      //         stringstream convert;
      // 	model.storage->ClearSample(level+2);  
      // 	convert.str(string());
      //         convert << model.parameter->run_id << "/" << model.parameter->run_id << "*." << level+2 << ".*";;
      //         string remove_file = model.parameter->storage_dir + convert.str();
      // 	remove(remove_file.c_str()); 
      // }

      //if (model.K(model.energy_level) < 1.0000001) break;
    }

  delete []sPackage; 
  delete []rPackage; 
}
