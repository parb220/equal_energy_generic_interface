#include <mpi.h>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include "dw_dense_matrix.hpp"
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CSampleIDWeight.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"
#include "prcsn.h"

#include "slave_computing.h"

using namespace std;

void slave_computing(CEquiEnergyModel &model, const CSampleIDWeight &mode)
{
	int my_rank, nNode; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &nNode);  
	MPI_Status status; 
	
	double *rPackage = new double [N_MESSAGE], *sPackage = new double [N_MESSAGE];    
	int group_index; 

	while (1)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
		if (status.MPI_TAG == END_TAG)
		{
		  // cout << "* END_TAG\n";
			delete []rPackage; 
			delete []sPackage; 
			exit(0); 
		}
		else if (status.MPI_TAG == BINNING_INFO)
		{
		  // cout << "* BINNING_INFO\n";
			model.energy_stage = (int)(rPackage[LEVEL_INDEX]); 
			int number_ring = (int)(rPackage[RESERVE_INDEX]); 
			model.storage-> ResizeBin(model.energy_stage, number_ring); 
			for (int i=0; i<number_ring; i++)
				model.storage->SetEnergyLowerBound(model.energy_stage, i, rPackage[RESERVE_INDEX+i+1]);
			if (model.energy_stage+1 <= model.parameter->number_energy_stage)
				model.storage->ClearBin(model.energy_stage+1); 
			sPackage[RETURN_INDEX_1] = 0;
                        sPackage[RETURN_INDEX_2] = 0;
		}		
		else if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST || status.MPI_TAG == SIMULATION_TAG || status.MPI_TAG == SIMULATION_PRIOR_TAG || status.MPI_TAG == SCALE_MATRIX_FIT_TAG) 
		{	
		  // cout << "* TUNE_TAG_SIMULATION_FIRST || SIMULATION_TAG || SIMULATION_PRIOR_TAG\n";
			model.energy_stage = (int)(rPackage[LEVEL_INDEX]);
			group_index = (int)(rPackage[GROUP_INDEX]); 
			int nGroup = (int)(rPackage[GROUP_NUMBER_INDEX]); 
			// model.timer_when_started = group_index; 
		
			if (!GetCommunicationParameter(rPackage, N_MESSAGE, model.parameter))
			{
				cout << "GetCommunicationParameter() : Error occurred.\n"; 
				abort(); 
			}
			model.lambda = model.parameter->lambda[model.energy_stage];
			std::vector<int> nJump = ExecutingSimulationTask(model, my_rank, group_index, nGroup,  mode, status.MPI_TAG); 
			sPackage[RETURN_INDEX_1] = nJump[0];
                        sPackage[RETURN_INDEX_2] = nJump[1];
		}
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
	}
}

