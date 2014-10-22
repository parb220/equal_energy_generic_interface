#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <time.h>
#include "CEquiEnergyModel.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

using namespace std;

void master_deploying(int nNode, CEquiEnergyModel &model)
{	
  time_t start_time=time((time_t*)NULL);

  if ((model.parameter->N > 0) && (model.parameter->Gn > 0))
    {
      // Write parameters
      model.WriteParameters();

      // Initial distribution
      if (model.parameter->first_stage == 0)
	DispatchInitialSimulation(nNode, model);

      // Simulation
      DispatchTuneSimulation(nNode, model, false); 
    }

  cout << "Ellapsed time: " << difftime(time((time_t*)NULL),start_time) << " seconds " << endl;

  // tell all the slaves to exit by sending an empty messag with 0 simulation length 
  double *sMessage= new double [N_MESSAGE];  
  for (int i=1; i<nNode; i++)
    MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
  delete [] sMessage; 
}

void DispatchInitialSimulation(int nNode, CEquiEnergyModel &model)
{
  MPI_Status status;
  double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];

  // diagnostic file
  fstream diagnostic_file;
  model.OpenFile(diagnostic_file,true,"Diagnostic",0);

  TDenseVector mean=Ones(model.nParameters);
  TDenseMatrix variance=Identity(model.nParameters);

  // send generate draws from initial distribution command to compute nodes
  model.SetupMasterNode(mean,variance);
  for (int i=1; i<nNode; i++)
    {
      sPackage[STAGE_INDEX] = 0.0;
      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG_INITIAL_INDEPENDENT, MPI_COMM_WORLD);		
    }
  for (int i=1; i<nNode; i++)
    {
      MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG_INITIAL_INDEPENDENT, MPI_COMM_WORLD, &status);		
    }
  model.storage->FinalizeMasterNode(0);

  bool not_done=true;
  int iterations=0;
  do
    {
     // send TUNE_TAG_INITIAL command to compute nodes
      model.SetupMasterNode(1);
      double scale = 0.0;
      for (int i=1; i<nNode; i++)
	{
	  sPackage[STAGE_INDEX] = 1.0; 
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG_INITIAL, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_INITIAL, MPI_COMM_WORLD, &status);	
	  scale+=rPackage[SCALE_INDEX];
	}
      scale/=(double)(nNode - 1);

      // send SIMULATION_TAG_INITIAL_METROPOLIS command to compute nodes
      for (int i=1; i<nNode; i++)
	{
	  sPackage[STAGE_INDEX] = 1.0;
	  sPackage[SCALE_INDEX] = scale;
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG_INITIAL_METROPOLIS, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG_INITIAL_METROPOLIS, MPI_COMM_WORLD, &status);
	}
      model.storage->FinalizeMasterNode(1);

      // send SIMULATION_TAG_INITIAL_INDEPENDENT command to compute nodes
      model.ComputeUnweightedMeanVariance(1,mean,variance);
      model.SetupMasterNode(mean,variance);
      for (int i=1; i<nNode; i++)
	{
	  sPackage[STAGE_INDEX] = 0.0; 
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG_INITIAL_INDEPENDENT, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG_INITIAL_INDEPENDENT, MPI_COMM_WORLD, &status);		
	}
      model.storage->FinalizeMasterNode(0);

      // write to diagnostic file
      model.SetupMasterNode(1);
      diagnostic_file << "Effective sample size = " << model.ESS(1) << endl;
      diagnostic_file << "Log of the integral of the kernel = " << model.LogIntegralKernel(1) << endl;
      diagnostic_file << "Scale = " << scale << endl;
      diagnostic_file << "Sqrt Diagonal =" << endl << model.InitialSqrtDiagonal << endl;
      diagnostic_file << "Orthonormal Directions =" << endl << model.InitialOrthonormalDirections << endl;
      diagnostic_file << "Mean =" << endl << model.InitialCenter << endl;  
      diagnostic_file << "/////////////////////////////////////////////////////////////////////////////////" << endl << endl;

      // exit if ESS is large enough or too many iterations have been performed
      iterations++;
      if ((iterations > 1) && ((model.ESS(1) >= model.parameter->min_ess) || (iterations > 5)))
      	{
      	  not_done=false;
      	}
    }
  while (not_done);

  diagnostic_file.close();

  delete[]sPackage;
  delete[]rPackage;
}

void DispatchTuneSimulation(int nNode, CEquiEnergyModel &model, bool save_space_flag)
{
  double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
  MPI_Status status; 

  int stage=(model.parameter->first_stage == 0) ? 1 : model.parameter->first_stage;
  for (double lambda=0.0; lambda < 0.9999999; lambda=model.lambda(stage++))
    {
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Send out tuning jobs
      cout << "Dispatching tuning tasks to nodes - stage " << stage << endl;
      model.SetupMasterNode(stage);
      double scale=0.0;
      for (int i=1; i<nNode; i++)
	{
	  sPackage[STAGE_INDEX]=(double)stage;
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG, MPI_COMM_WORLD, &status);	
	  scale+=rPackage[SCALE_INDEX];
	}
      model.scale=scale/=(double)(nNode - 1);

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Send out simulation jobs
      cout << "Dispatching simulation tasks to nodes - stage " << stage << endl; 	
      for (int i=1; i<nNode; i++)
	{
	  sPackage[STAGE_INDEX]=(double)stage;
	  sPackage[SCALE_INDEX] = scale;
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG, MPI_COMM_WORLD);		
	}	
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG, MPI_COMM_WORLD, &status);
	}
      model.storage->FinalizeMasterNode(stage);
      model.WriteDiagnosticMasterNode(stage);

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
    }

  model.storage->WriteASCII_draws(stage-1);

  delete []sPackage; 
  delete []rPackage; 
}
