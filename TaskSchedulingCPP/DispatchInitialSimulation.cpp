

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
  MPI_Status status;
  double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];
  sPackage[LEVEL_INDEX] = (double)(model.energy_level);

  // diagnostic file
  fstream diagnostic_file;
  model.OpenFile(diagnostic_file,true,"Diagnostic",model.energy_level);

  // create first initial distribution file
  model.energy_level=model.parameter->number_energy_level;
  model.InitialOrthonormalDirections.Identity(model.nParameters);
  model.InitialSqrtDiagonal.Ones(model.nParameters);
  model.InitialCenter.Zeros(model.nParameters);
  model.WriteInitialDistribution();

  // send generate draws from initial distribution command to compute nodes
  for (int i=1; i<nNode; i++)
    {
      sPackage[GROUP_INDEX] = (double)(i-1); 
      MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG_INITIAL_INDEPENDENT, MPI_COMM_WORLD);		
    }
  for (int i=1; i<nNode; i++)
    {
      MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG_INITIAL_INDEPENDENT, MPI_COMM_WORLD, &status);		
    }
  model.storage->ConsolidateDraws(model.energy_level, model.energy_level);

  bool not_done=true;
  int iterations=0;
  do
    {
      double scale = 0.0;
     // send tune command to compute nodes
      for (int i=1; i<nNode; i++)
	{
	  sPackage[GROUP_INDEX] = (double)(i-1); 
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG_INITIAL, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_INITIAL, MPI_COMM_WORLD, &status);	
	  scale+=rPackage[SCALE_INDEX];
	}
      scale/=(double)(nNode - 1);

      // send metropolis simulate command to compute nodes
      for (int i=1; i<nNode; i++)
	{
	  sPackage[GROUP_INDEX] = (double)(i-1); 
	  sPackage[SCALE_INDEX] = scale;
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG_INITIAL_METROPOLIS, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG_INITIAL_METROPOLIS, MPI_COMM_WORLD, &status);
	}
      cout << "consolidating files\n";
      model.storage->ConsolidateDraws(model.energy_level-1, model.energy_level-1);
      model.storage->SetupForInput(model.energy_level-1,1,model.parameter->max_energy);
      model.storage->ComputeVarianceDecomposition(model.InitialOrthonormalDirections,model.InitialSqrtDiagonal);
      model.storage->ComputeMean(model.InitialCenter);
      model.WriteInitialDistribution();

      // send generate draws from initial distribution command to compute nodes
      for (int i=1; i<nNode; i++)
	{
	  sPackage[GROUP_INDEX] = (double)(i-1); 
	  MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SIMULATION_TAG_INITIAL_INDEPENDENT, MPI_COMM_WORLD);		
	}
      for (int i=1; i<nNode; i++)
	{
	  MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SIMULATION_TAG_INITIAL_INDEPENDENT, MPI_COMM_WORLD, &status);		
	}
      model.storage->ConsolidateDraws(model.energy_level, model.energy_level);

      // write to diagnostic file
      model.storage->SetupForInput(model.energy_level,1,model.parameter->max_energy);
      double ess = model.storage->ComputeEffectiveSampleSize();
      diagnostic_file << "Effective sample size = " << ess << endl;
      diagnostic_file << "Log of the integral of the kernel = " << model.storage->ComputeLogIntegralKernel() << endl;
      diagnostic_file << "Scale = " << scale << endl;
      diagnostic_file << "Sqrt Diagonal =" << endl << model.InitialSqrtDiagonal << endl;
      diagnostic_file << "Orthonormal Directions =" << endl << model.InitialOrthonormalDirections << endl;
      diagnostic_file << "Mean =" << endl << model.InitialCenter << endl;  
      diagnostic_file << "/////////////////////////////////////////////////////////////////////////////////" << endl << endl;

      // exit if ESS is large enough or too many iterations have been performed
      iterations++;
      if ((iterations > 1) && ((ess >= model.parameter->min_ess) || (iterations > 5)))
	{
	  not_done=false;
	}
    }
  while (not_done);

  diagnostic_file.close();

  delete[]sPackage;
  delete[]rPackage;
}

