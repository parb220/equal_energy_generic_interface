#include <vector>
#include <sstream>
#include <glob.h>
#include <mpi.h> 
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
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


double DispatchSimulation(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, int level, int message_tag)
{
	double max_log_posterior =-1.0e300, received_log_posterior; 
	double *sPackage = new double [N_MESSAGE]; 
	// burn_in_length: 0.1*simulation_length_per_node or 5000, whichever is larger
	sPackage[FREQ_INDEX] = parameter.deposit_frequency; 
       	sPackage[LEVEL_INDEX] = level;
	sPackage[H0_INDEX] = parameter.h0; 

	size_t nNode=0; 
	for (int m=0; m<nodeGroup.size(); m++)
		nNode += nodeGroup[m].size(); 
	size_t simulation_length_per_node; 
	for (int i=0; i<nodeGroup.size(); i++)
	{
		if (message_tag == TUNE_TAG_SIMULATION_FIRST)
			simulation_length_per_node = (size_t)ceil((double)simulation_length/(double)nodeGroup[i].size()); 
		else
			simulation_length_per_node = (size_t)ceil((double)simulation_length/(double)nNode); 
		for (int j = 0; j<nodeGroup[i].size(); j++)
		{
			sPackage[LENGTH_INDEX] = simulation_length_per_node; 
			sPackage[BURN_INDEX] = (simulation_length_per_node*parameter.deposit_frequency)/10 >= 5000 ? (simulation_length_per_node*parameter.deposit_frequency)/10 : 5000; 
			sPackage[GROUP_INDEX] = i; 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodeGroup[i][j], message_tag, MPI_COMM_WORLD);
		}
	}
	delete [] sPackage;

	MPI_Status status;
	double *rPackage = new double [N_MESSAGE];
	for (int i=0; i<nodeGroup.size(); i++)
	{
		for (int j=0; j<nodeGroup[i].size(); j++)
		{
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, message_tag, MPI_COMM_WORLD, &status); 
			received_log_posterior = rPackage[H0_INDEX]; 
			max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior; 
		}
	}
	delete [] rPackage;

	// Consolidate partial storage files
	storage.consolidate(level); 

	// Consolidate variance file
	if (message_tag == TUNE_TAG_SIMULATION_FIRST || message_tag == TUNE_TAG_SIMULATION_SECOND)
	{
		stringstream convert; 
		string input_file_pattern, output_file; 
		if (message_tag == TUNE_TAG_SIMULATION_FIRST)
		{
			for (int i=0; i<nodeGroup.size(); i++)
			{
				convert.str(string()); 
				convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << level << "." << i << ".*";
				input_file_pattern = parameter.storage_dir + convert.str();
				
				convert.str(string());
                        	convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << level << "." << i;
				output_file = parameter.storage_dir + convert.str(); 
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
                        convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << level << ".*.*";
			input_file_pattern = parameter.storage_dir + convert.str();
			convert.str(string());
			convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << level;
			output_file = parameter.storage_dir + convert.str(); 
			if (!ConsolidateSampleForCovarianceEstimation(input_file_pattern, output_file))
                       	{
                       		cerr << "ConsolidateSampleForCovarianceEstimation() : Error.\n";
                       		abort();
                       	}
		}
	}
	return max_log_posterior; 
}
