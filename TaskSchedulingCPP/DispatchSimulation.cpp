#include <vector>
#include <sstream>
#include <glob.h>
#include <mpi.h> 
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "dw_math.h"
#include "CEquiEnergyModel.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "storage_parameter.h"
#include "mpi_parameter.h"

using namespace std;  

size_t glob(vector<string> &filename, const string &pattern)
{
        glob_t glob_result;
        glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
        if (glob_result.gl_pathc > 0)
        {
		filename.resize(glob_result.gl_pathc); 
                for (int i=0; i<glob_result.gl_pathc; i++)
                        filename[i] = string(glob_result.gl_pathv[i]);
	}
        globfree(&glob_result);
        return filename.size();
}

bool ConsolidateSampleForCovarianceEstimation(const string &file_pattern, const string &variance_file)
{
        vector <string> filenames_merge; 
	size_t number_file = glob(filenames_merge, file_pattern); 
	if (number_file == 0)
		return true; 
        ofstream out_file(variance_file.c_str(), ios::out | ios::binary);
        if (!out_file)
                return false;
        ifstream input_file;
        int fail_counter =0;
        for (int i=0; i<filenames_merge.size(); i++)
        {
                input_file.open(filenames_merge[i].c_str(), ios::in | ios::binary);
                if (!input_file)
                {
                        fail_counter ++;
                        cerr << "ConsolidateSampleForCovarianceEstimation(): Error in opening " << filenames_merge[i] << "for reading.\n";
                        continue;
                }
                out_file << input_file.rdbuf();
                out_file.flush();
                input_file.close();
                remove(filenames_merge[i].c_str());
        }
        out_file.close();
        if (fail_counter < filenames_merge.size() ) 
                return true;
	else 
		return false; 
}


void DispatchSimulation(int nNode, int nInitial, CEquiEnergyModel &model, size_t simulation_length, int level, int message_tag)
{
	double *sPackage = new double [N_MESSAGE]; 
	double *rPackage = new double [N_MESSAGE]; 
	// burn_in_length: 0.1*simulation_length_per_node or 5000, whichever is larger
	sPackage[thin_INDEX] = model.parameter->thin; 
	sPackage[THIN_INDEX] = model.parameter->THIN;
       	sPackage[LEVEL_INDEX] = level;
	sPackage[PEE_INDEX] = model.parameter->pee; 

	size_t simulation_length_per_node; 
	if (message_tag == TUNE_TAG_SIMULATION_FIRST)
		simulation_length_per_node = (size_t)ceil((double)simulation_length/(double)(nNode-1)) > 1000 ? (size_t)ceil((double)simulation_length/(double)(nNode-1)) : 1000; 
	else 
		simulation_length_per_node = (size_t)ceil((double)simulation_length/(double)(nInitial*(nNode-1))) > 1000 ? (size_t)ceil((double)simulation_length/(double)(nInitial*(nNode-1))) : 1000; 
	MPI_Status status;
	
	vector<int> available_node(nNode-1); 
	for (int i=0; i<(int)(available_node.size()); i++)
		available_node[i] = i+1; 
	
	int iInitial =0, cumulative_length = 0; 
	while (iInitial < nInitial)
	{
		while (cumulative_length < simulation_length)
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
	}
	
	for (int j=0; j<(nNode-1)-(int)available_node.size(); j++)
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, message_tag, MPI_COMM_WORLD, &status);

	/*for (int iInitial=0; iInitial<nInitial; iInitial++)
	{


		for (int j = 1; j<nNode; j++)
		{
			sPackage[LENGTH_INDEX] = simulation_length_per_node; 
			sPackage[BURN_INDEX] = model.parameter->burn_in_length; //(simulation_length_per_node*model.parameter->thin)/10 >= BURN_IN_LENGTH ? (simulation_length_per_node*model.parameter->thin)/10 : BURN_IN_LENGTH; 
			sPackage[GROUP_INDEX] = iInitial; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, j, message_tag, MPI_COMM_WORLD);
		}

		for (int j=1; j<nNode; j++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, message_tag, MPI_COMM_WORLD, &status); 
	}*/
	delete [] sPackage;
	delete [] rPackage;

	// Consolidate partial storage files
	model.storage->ClearStatus(level); 
	model.storage->consolidate(level); 

	// Consolidate variance file
	if (message_tag == TUNE_TAG_SIMULATION_FIRST || message_tag == TUNE_TAG_SIMULATION_SECOND)
	{
		stringstream convert; 
		string input_file_pattern, output_file; 
		if (message_tag == TUNE_TAG_SIMULATION_FIRST)
		{
			for (int iInitial=0; iInitial<nInitial; iInitial++)
			{
				convert.str(string()); 
				convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << level << "." << iInitial << ".*";
				input_file_pattern = model.parameter->storage_dir + convert.str();
				
				convert.str(string());
                        	convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << level << "." << iInitial;
				output_file = model.parameter->storage_dir + convert.str(); 
				if (!ConsolidateSampleForCovarianceEstimation(input_file_pattern, output_file))
                        	{
                                	cerr << "ConsolidateSampleForCovarianceEstimation() : Error.\n";
                                	abort();
                        	}
			}
		}
		else
		{
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << level << ".*.*";
			input_file_pattern = model.parameter->storage_dir + convert.str();
			convert.str(string());
			convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << level;
			output_file = model.parameter->storage_dir + convert.str(); 
			if (!ConsolidateSampleForCovarianceEstimation(input_file_pattern, output_file))
                       	{
                       		cerr << "ConsolidateSampleForCovarianceEstimation() : Error.\n";
                       		abort();
                       	}
		}
	}
}
