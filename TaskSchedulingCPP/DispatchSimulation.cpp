#include <vector>
#include <sstream>
#include <mpi.h> 
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "dw_math.h"
#include "CSampleIDWeight.h"
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "storage_parameter.h"
#include "mpi_parameter.h"

using namespace std;  

vector<string> glob(const string &pattern); 

void DispatchSimulation_PriorProbability(int nNode, int simulation_length)
{
	double *sPackage = new double [N_MESSAGE];
        double *rPackage = new double [N_MESSAGE];

	int simulation_length_per_node =  (int)ceil((double)simulation_length/(double)(nNode-1)); 
	int total_length = 0; 
	for (int iNode=1; iNode<nNode; iNode++)
	{
		sPackage[LENGTH_INDEX] = simulation_length_per_node;
		total_length += sPackage[LENGTH_INDEX]; 
		sPackage[GROUP_INDEX] = iNode; 
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, iNode, PRIOR_PROB_TAG, MPI_COMM_WORLD);
	}
	int accpt_length = 0; 
	MPI_Status status;
	for (int iNode=1; iNode<nNode; iNode++)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, PRIOR_PROB_TAG, MPI_COMM_WORLD, &status);
		accpt_length += rPackage[LENGTH_INDEX]; 
	}
	cout << "log prior probability constant " << setprecision(20) << log((double)accpt_length/(double)total_length) << endl; 
}

std::vector<CSampleIDWeight> DispatchSimulation(int nNode, int nInitial, CEquiEnergyModel &model, int simulation_length, int stage, int message_tag)
{
	double *sPackage = new double [N_MESSAGE]; 
	double *rPackage = new double [N_MESSAGE]; 
	sPackage[thin_INDEX] = model.parameter->thin; 
	sPackage[THIN_INDEX] = model.parameter->THIN;
       	sPackage[LEVEL_INDEX] = stage;
	sPackage[PEE_INDEX] = model.parameter->pee; 
	std::vector<CSampleIDWeight> samples(0); 

	MPI_Status status;
	
	vector<int> available_node(nNode-1); 
	for (int i=0; i<(int)(available_node.size()); i++)
		available_node[i] = i+1; 
	
	std::vector<int> length_task;
	std::vector<int> group_index_task;
	int length_per_initial = (int)ceil((double)simulation_length/(double)nInitial);
	if (length_per_initial > TASK_LENGTH)
	{
		for (int iInitial=0; iInitial<nInitial; iInitial++)
		{
			for (int k=0; k<length_per_initial/TASK_LENGTH; k++)
			{
				length_task.push_back(TASK_LENGTH);
				group_index_task.push_back(iInitial);
			}
			if (length_per_initial%TASK_LENGTH)
			{
				length_task.push_back(length_per_initial % TASK_LENGTH);
				group_index_task.push_back(iInitial);
			}
		}

	}
	else
	{
		for (int iInitial=0; iInitial<nInitial; iInitial++)
		{
			length_task.push_back(length_per_initial);
			group_index_task.push_back(iInitial);
		}	
	}

	for (int iTask = 0; iTask<(int)(length_task.size()); iTask ++)
	{
		sPackage[LENGTH_INDEX] = length_task[iTask];  
		sPackage[BURN_INDEX] = model.parameter->burn_in_length;
		sPackage[GROUP_INDEX] = group_index_task[iTask]; 
		if (available_node.empty())
		{
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, message_tag, MPI_COMM_WORLD, &status); 
			available_node.push_back(status.MPI_SOURCE); 
		}
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, available_node.back(), message_tag, MPI_COMM_WORLD);
		available_node.pop_back(); 
	}

	/*int simulation_length_per_node, simulation_length_check; 
	// simulation_length_per_node = (int)ceil((double)simulation_length/(double)(nInitial*(nNode-1))) > 1000 ? (int)ceil((double)simulation_length/(double)(nInitial*(nNode-1))) : 1000; 
	simulation_length_per_node = (int)ceil((double)simulation_length/(double)(nInitial*(nNode-1))) > 1 ?  (int)ceil((double)simulation_length/(double)(nInitial*(nNode-1))): 1 ; 
	simulation_length_check = (int)ceil((double)simulation_length/(double)(nInitial)); 

	
	int iInitial =0, cumulative_length = 0; 
	while (iInitial < nInitial)
	{
		while (cumulative_length < simulation_length_check)
		{
			sPackage[LENGTH_INDEX] = simulation_length_per_node;
			sPackage[BURN_INDEX] = model.parameter->burn_in_length;
			sPackage[GROUP_INDEX] = iInitial; 
			if (available_node.empty())
			{
				MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, message_tag, MPI_COMM_WORLD, &status); 
				available_node.push_back(status.MPI_SOURCE); 
			}
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, available_node.back(), message_tag, MPI_COMM_WORLD);
			available_node.pop_back(); 
			cumulative_length += simulation_length_per_node; 	
		}
		cumulative_length = 0; 
		iInitial ++; 
	}*/
	
	for (int j=0; j<(nNode-1)-(int)available_node.size(); j++)
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, message_tag, MPI_COMM_WORLD, &status);

	// binning
	if (message_tag != TUNE_TAG_SIMULATION_FIRST) 
	{
		model.storage->ClearStatus(stage);	
		model.storage->consolidate(stage); 
		samples = model.storage->binning_equal_size(stage, model.parameter->number_striation);
		model.storage->finalize(stage);
        	model.storage->ClearDepositDrawHistory(stage);

		sPackage[LEVEL_INDEX] = (double)stage;
        	sPackage[RESERVE_INDEX] = (double)(model.storage->GetNumber_Bin(stage));
        	for (int i=0; i<model.storage->GetNumber_Bin(stage); i++)
        		sPackage[RESERVE_INDEX + i + 1] = model.storage->GetEnergyLowerBound(stage, i);

        	for (int i=1; i<nNode; i++)
        		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, BINNING_INFO, MPI_COMM_WORLD);

        	for (int i=1; i<nNode; i++)
        		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, BINNING_INFO, MPI_COMM_WORLD, &status);
	}
	else // Consolidate variance file
	{
		stringstream convert; 
		string input_file_pattern, output_file_name;
		convert.str(string());
                convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << stage << ".*.*";
		input_file_pattern = model.parameter->storage_dir + convert.str();

		convert.str(string());
		convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << stage; 
		output_file_name = model.parameter->storage_dir + convert.str();

		vector<string> input_file_name = glob(input_file_pattern);  

		if (!model.ConsolidateSampleForCovarianceEstimation(input_file_name, output_file_name))
                {
                	cerr << "ConsolidateSampleForCovarianceEstimation() : Error.\n";
                	abort();
                }
	}
	
	delete [] sPackage;
	delete [] rPackage;
	return samples;
}
