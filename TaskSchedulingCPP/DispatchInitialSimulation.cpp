

#include <vector>
#include <cmath>
#include <mpi.h>
#include <glob.h>
#include <sstream>
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

using namespace std; 

void DispatchInitialSimulation(int nNode, CEquiEnergyModel &model)
{
  int nParameters=model.nParameters;
  model.energy_level=model.parameter->number_energy_level;
  model.node=-1;

  MPI_Status status;
  double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];
  sPackage[LEVEL_INDEX] = (double)(model.energy_level);

  // diagnostic file
  fstream diagnostic_file;
  model.OpenFile(diagnostic_file,"Diagnostic",model.energy_level,true);

  // create initialization file
  model.SetupFromPreviousLevel(model.energy_level);
 
  bool not_done=true;
  int iterations=0;
  do
    {
      // send out tune command to compute nodes
      for (int i=1; i<nNode; i++)
	{
	  sPackage[GROUP_INDEX] = i-1; 
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG_INITIAL, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_INITIAL, MPI_COMM_WORLD, &status);		
	}

      // consolidate scales and create initialization file
      model.scale = model.ConsolidateScales(model.energy_level);
      model.WriteInitializationFile();
      diagnostic_file << "Temperature = " << model.parameter->max_energy << endl;
      diagnostic_file << "Scale = " << model.scale << endl;

      // send out simulate command to compute nodes
      for (int i=1; i<nNode; i++)
	{
	  sPackage[GROUP_INDEX] = i-1; 
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG_INITIAL, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG_INITIAL, MPI_COMM_WORLD, &status);
	}

      // read draws, set , compute ESS
      double ESS = model.AnalyzeInitialDraws();
      diagnostic_file << "ESS = " << ESS << endl;

      cout << "ESS (" << model.parameter->min_ess << ") = " << ESS << endl << endl;

      // exit if ESS is large enough or too many iterations have been performed
      not_done=false;
      iterations++;
      if ((iterations > 1) && ((ESS >= model.parameter->min_ess) || (iterations > 5)))
	{
	  not_done=false;
	}
    }
  while (not_done);

  // simulate and save top level draws
  TDenseVector x;
  CSampleIDWeight y;
  int length=model.parameter->simulation_length;
  TDenseVector log_kernel_target(length), log_density_proposal(length);
  for (int ii=0; ii < length; ii++)
    {
      log_density_proposal(ii)=model.InitialDraw(x);
      y.data=x;
      y.weight=model.target->LogPosterior(x);
      y.id=(int)(time(NULL)-model.timer_when_started);
      y.calculated=true;
      model.SaveSampleToStorage(y); 
      log_kernel_target(ii)=y.weight/model.parameter->max_energy;
    }

  diagnostic_file << "ESS = " << EffectiveSampleSize(log_kernel_target,log_density_proposal) << endl;
  diagnostic_file << "MDD = " << ComputeLogKernelIntegral(log_kernel_target,log_density_proposal) << endl;
  diagnostic_file.close();

  // clean up
  model.storage->ClearStatus(model.energy_level);
  model.storage->binning_equal_size(model.energy_level, model.parameter->number_rings);
  model.storage->finalize(model.energy_level);
  model.storage->ClearDepositDrawHistory(model.energy_level); 

  sPackage[LEVEL_INDEX] = (double)(model.energy_level);
  sPackage[RESERVE_INDEX] = (double)(model.storage->GetNumber_Bin(model.energy_level)); 
  for (int i=0; i < model.storage->GetNumber_Bin(model.energy_level); i++)
    sPackage[RESERVE_INDEX + i + 1] = model.storage->GetEnergyLowerBound(model.energy_level, i); 
		
  for (int i=1; i<nNode; i++)
    MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, BINNING_INFO, MPI_COMM_WORLD); 

  for (int i=1; i<nNode; i++)
    MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, BINNING_INFO, MPI_COMM_WORLD, &status);

  delete[]sPackage;
  delete[]rPackage;
}

