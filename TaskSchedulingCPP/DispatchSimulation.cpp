#include <sstream>
#include <glob.h>
#include "master_deploying.h"
#include "storage_parameter.h"

using namespace std;  

size_t glob(vector<string> &filename, const string &pattern)
{
        glob_t glob_result;
        glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
        if (glob_result.gl_pathc > 0)
        {
		filename.resize(glob_result.gl_pathc); 
                for (unsigned int i=0; i<glob_result.gl_pathc; i++)
                        filename[i] = string(glob_result.gl_pathv[i]);
	}
        globfree(&glob_result);
        return filename.size();
}

bool ConsolidateSampleForCovarianceEstimation(const CEESParameter &parameter, unsigned int level)
{
        stringstream convert;
        convert.str(string());
        convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << level << ".*";
        string file_pattern = parameter.storage_dir + convert.str();

        vector <string> filenames_merge; 
        convert.str(string());
	size_t number_file = glob(filenames_merge, file_pattern); 
	if (number_file == 0)
		return true; 
        convert << parameter.run_id << "/" << parameter.run_id << VARIANCE_SAMPLE_FILE_TAG << level;
        string variance_file = parameter.storage_dir + convert.str();
        ofstream out_file(variance_file.c_str(), ios::out | ios::binary);
        if (!out_file)
                return false;
        ifstream input_file;
        unsigned int fail_counter =0;
        for (unsigned i=0; i<filenames_merge.size(); i++)
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
}


double DispatchSimulation(const vector<unsigned int> &nodePool, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, unsigned int level, int message_tag)
{
	double *sPackage = new double [N_MESSAGE]; 
	size_t simulation_length_per_node = (size_t)((double)simulation_length/(double)nodePool.size()); 
	sPackage[LENGTH_INDEX] = simulation_length_per_node; 
	// burn_in_length: 0.1*simulation_length_per_node or 5000, whichever is larger
	sPackage[BURN_INDEX] = (simulation_length_per_node*parameter.deposit_frequency)/10 >= 5000 ? (simulation_length_per_node*parameter.deposit_frequency)/10 : 5000; 
	sPackage[FREQ_INDEX] = parameter.deposit_frequency; 
       	sPackage[LEVEL_INDEX] = level;
	sPackage[H0_INDEX] = parameter.h0; 

	for (unsigned int i=0; i<nodePool.size(); i++)
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, nodePool[i], message_tag, MPI_COMM_WORLD);
	delete [] sPackage;

	MPI_Status status;
	double *rPackage = new double [N_MESSAGE];
	double max_log_posterior =-1.0e300, received_log_posterior; 
	for (unsigned int i=0; i<nodePool.size(); i++)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, message_tag, MPI_COMM_WORLD, &status); 
		received_log_posterior = rPackage[H0_INDEX]; 
		max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior; 
	}
	delete [] rPackage;

	// Consolidate partial storage files
	int start_bin = parameter.BinIndex_Start(level); 
	int end_bin = parameter.BinIndex_End(level); 
	storage.consolidate(start_bin, end_bin); 

	// Consolidate variance file
	if (message_tag == TUNE_TAG_SIMULATION_FIRST || message_tag == TUNE_TAG_SIMULATION_SECOND)
	{
		if (!ConsolidateSampleForCovarianceEstimation(parameter, level))
		{
			cerr << "ConsolidateSampleForCovarianceEstimation() : Error.\n";
                	abort();
		}
	}
	return max_log_posterior; 
}
