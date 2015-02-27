#include <cmath>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "dw_rand.h"
#include "dw_math.h"
#include "CSampleIDWeight.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"
#include "storage_parameter.h"
#include "mpi_parameter.h"
#include "dw_matrix.h"
#include "mdd.hpp"
#include "mdd_function.h"
#include "master_deploying.h"
#include "EstimateLogMDD.hpp"

using namespace std; 

vector<string> glob(const string &pattern);

void DispatchTuneSimulation(int nNode, int nInitial, CEquiEnergyModel &model,const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag)
{
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
	MPI_Status status; 

	size_t estimation_length; 

	// diagnostic statistics
	// logMDD[i][0], constant of integration using ellipticals as proposals
	// logMDD[i][1], constant of integration using draws from the previous stage and importance.
	// logMDD[i][2], constant of integration using draws from the previous stage as the proposals, and logMDD[i-1][0] as the estimate for the previous stage 
	// logMDD[i][3], constant of integration using draws from the previous stage as the proposals, and logMDD[i-1][3] as the estimate for the previous stage
	//
	vector<vector<double> > logMDD(model.parameter->number_energy_stage+1, vector<double>(4, 0.0)); 
	/* vector<double> consistency(model.parameter->number_energy_stage, 0.0);
	vector<double> average_consistency(model.parameter->number_energy_stage, 0.0); 
	vector<double> std_consistency(model.parameter->number_energy_stage, 0.0);
	vector<double> LB_ESS(model.parameter->number_energy_stage, 0.0); */

	vector<CSampleIDWeight> posterior, proposal;
	int data_size = model.current_sample.GetSize_Data();
	bool unstructured = true;
	
	for (int stage=model.parameter->highest_stage; stage>=model.parameter->lowest_stage; stage--)
	{
		/////////////////////////////////////////////////////////////////////////////////
		// Highest +1 stage 
		if (stage == model.parameter->number_energy_stage-1)
		{
			// HighestPlus1Stage(nNode, nInitial, model);
			HighestPlus1Stage_Prior(nNode, nInitial, model); // Sample from prior
			posterior.clear(); 
			if (!model.storage->DrawAllSample(stage+1, posterior, unstructured, data_size))
                        {
                                cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
                                abort();
                        }

			// logMDD[stage+1][1]: logMDD using elliptical for draws from prior distribution
                        logMDD[stage+1][0] = logMDD[stage+1][1] = logMDD[stage+1][2] = logMDD[stage+1][3] = LogMDD(posterior, model, model.parameter->t[stage+1], USE_TRUNCATED_POWER, PRIOR_ONLY);
			cout << setprecision(20) << "logMDD at stage " << stage+1 << ": " << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << "\t" << logMDD[stage+1][2] << "\t" << logMDD[stage+1][3] << endl; 
		}

		////////////////////////////////////////////////////////////////////////////////
		// Starting points
		vector<CSampleIDWeight> start_points(nInitial); 
		if (model.storage->empty(stage+1) || !model.Initialize_WeightedSampling(nInitial, stage+1, start_points))
		{
			for (int i=0; i<nInitial; i++)
				start_points[i] = mode; 
		} 
		model.storage->ClearDepositDrawHistory(stage+1); 
		model.storage->ClearStatus(stage+1); 
		model.storage->RestoreForFetch(stage+1); 

		stringstream convert; 
		string start_point_file; 
        	ofstream output_file;
		for (int i=0; i<nInitial; i++)
		{
			convert.str(string());
        		convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << stage << "." << i;
        		start_point_file = model.parameter->storage_dir + convert.str();
        		output_file.open(start_point_file.c_str(), ios::binary|ios::out);
        		if (!output_file)
			{
                		cerr << "Error in writing to " << start_point_file << endl; 
				abort(); 	
			}
        		else
               			write(output_file, &(start_points[i]));
        		output_file.close();
		}

		/////////////////////////////////////////////////////////////////////////////////
		// Send out tuning jobs
		// Only need to tune on each slave node noce, because the results will be aggregated
		if (stage != model.parameter->number_energy_stage)
		{
			// Calculate B for AdaptiveAfterSimulation_WeightedSampling_OnePass
			
			// samples
			vector<CSampleIDWeight> samples;
                	if (!model.storage->DrawAllSample(stage+1, samples) || samples.empty() )
                	{
                        	cerr << "DispatchTuneSimulation() : Error in loading samples of previous stage.\n";
                        	abort();
                	}

			// weight
			vector<double> log_weight = model.Reweight(samples, stage, stage+1);
        		double log_weight_sum = log_weight[0];
        		for (int i=1; i<(int)log_weight.size(); i++)
                		log_weight_sum = AddLogs(log_weight_sum, log_weight[i]);
        		vector<double> weight(log_weight.size(), 0.0);
        		for (int i=0; i<(int)log_weight.size(); i++)
                		weight[i] = exp(log_weight[i] - log_weight_sum);
			
			// block_scheme
			convert.str(string());
        		convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_SCHEME;
        		string block_scheme_file_name = model.parameter->storage_dir  + convert.str();
			std::vector<TIndex> block_scheme = ReadBlockScheme(block_scheme_file_name);
			if (block_scheme.empty())
		                block_scheme.push_back(TIndex(0, samples[0].data.dim-1));
			
			vector<TDenseMatrix>B_matrix = GetBlockMatrix_WeightedSampling(samples, weight, block_scheme); 

			// Bmatrix file
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << BMATRIX; 
			string bmatrix_file_name = model.parameter->storage_dir  + convert.str();
			if (!WriteBMatrixFile(bmatrix_file_name, B_matrix))
			{
				cerr << "DispatchTuneSimulation() : Error in writing BMatrix file.\n"; 
				abort(); 
			}
		}

		sPackage[LEVEL_INDEX] = stage; 
		sPackage[thin_INDEX] = model.parameter->thin; 
		sPackage[THIN_INDEX] = model.parameter->THIN; 
		sPackage[PEE_INDEX] = model.parameter->pee/(model.parameter->THIN/model.parameter->thin); 
		for (int i=1; i<nNode; i++)
		{
			sPackage[GROUP_INDEX] = dw_uniform_int(nInitial); 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD); 
		}
		
		for (int i=1; i<nNode; i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);

	
		/////////////////////////////////////////////////////////////////////////////////
		// Aggregate the scales obtained from computing nodes
		convert.str(string()); 
		convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << stage << ".*";
		string block_file_pattern = model.parameter->storage_dir + convert.str();
	        vector<string> block_file = glob(block_file_pattern);
		
		convert.str(string());
        	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << stage;
       	 	string block_file_name = model.parameter->storage_dir + convert.str();
		if (!model.metropolis->AggregateBlocksAndRemoveFiles(block_file, block_file_name))
        	{
                	cerr << "Error in reading " << block_file_pattern << " or writing " << block_file_name << endl;
                	abort();
        	}

		if (stage != model.parameter->number_energy_stage)
		{
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << BMATRIX;
                        string bmatrix_file_name = model.parameter->storage_dir  + convert.str();
			remove(bmatrix_file_name.c_str()); 
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// simualtion
		cout << "Simulation at " << stage << " for " << model.parameter->simulation_length << endl; 
		DispatchSimulation(nNode, nInitial, model, model.parameter->simulation_length, stage, SIMULATION_TAG);
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// logMDD
		/*if (stage == model.parameter->number_energy_stage -10)
		{
			posterior.clear(); 
 			if (!model.storage->DrawAllSample(stage+1, posterior, unstructured, data_size))
 			{
 				cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
 				abort();
 			}
 
			// logMDD[stage+1][1]: logMDD using elliptical for draws from prior distribution
			logMDD[stage+1][0] = logMDD[stage+1][1] = logMDD[stage+1][2] = logMDD[stage+1][3] = LogMDD(posterior, model, model.parameter->t[stage+1], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
			cout << setprecision(20) << "logMDD at stage " << stage+1 << ": " << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << "\t" << logMDD[stage+1][2] << "\t" << logMDD[stage+1][3] << endl; 
		}
		if (stage <= model.parameter->number_energy_stage -10)
		{*/
			// logMDD[stage][1]: logMDD using importance sampling of the draws from previous stage
			logMDD[stage][1] = LogMDD_Importance(posterior, model, stage+1, stage, logMDD[stage+1][1], LIKELIHOOD_HEATED); 
			proposal = posterior; 
			posterior.clear(); 
        		if (!model.storage->DrawAllSample(stage, posterior, unstructured, data_size))
        		{	
               			cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
               			abort();
        		}
	     		logMDD[stage][0] = LogMDD(posterior, model, model.parameter->t[stage], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
			logMDD[stage][2] = LogMDD(proposal, posterior, model, stage+1, stage, logMDD[stage+1][0]); 
			logMDD[stage][3] = LogMDD(proposal, posterior, model, stage+1, stage, logMDD[stage+1][3]); 
		
			cout << setprecision(20) << "logMDD at stage " << stage << ": " << logMDD[stage][0] << "\t" << logMDD[stage][1] << "\t" << logMDD[stage][2] << "\t" << logMDD[stage][3] << endl; 
		// }

		/*if (stage == model.parameter->number_energy_stage-1 )
			consistency[stage-1] = CheckConvergency(model, stage-1, stage, logMDD[stage][stage], average_consistency[stage-1], std_consistency[stage-1], LB_ESS[stage-1]);
		else
			consistency[stage-1] = CheckConvergency(model, stage-1, stage, consistency[stage], average_consistency[stage-1], std_consistency[stage-1], LB_ESS[stage-1]); 
                cout << "Convergency Measure at Stage " << stage-1 << ": " << setprecision(20) << consistency[stage-1] << "\t" << average_consistency[stage-1] << "\t" << std_consistency[stage-1]<< "\t" << LB_ESS[stage-1] << endl;  
		if (stage == 0)
		{
			logMDD[stage][stage] = EstimateLogMDD(model, stage, USE_TRUNCATED_POWER);
			cout << setprecision(20) << "logMDD at stage " << stage << ": " << logMDD[stage][stage] << endl; 
		} */
		
		// to save space, remove stage+1 samples
		if (save_space_flag  && stage > 0 ) //&& stage+1 < model.parameter->number_energy_stage-1 )
		{
			model.storage->ClearSample(stage+1);  
			
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << "*." << stage+1 << ".*";;
                        string remove_file_pattern = model.parameter->storage_dir + convert.str();
			vector<string> remove_file = glob(remove_file_pattern); 
                	for (int i=0; i<(int)remove_file.size(); i++)
                        	remove(remove_file[i].c_str());
			
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << "*.sample." << stage+1;
			remove_file_pattern = model.parameter->storage_dir + convert.str();
			remove_file = glob(remove_file_pattern); 
                        for (int i=0; i<(int)remove_file.size(); i++)
                                remove(remove_file[i].c_str());
		}

	}

	delete []sPackage; 
	delete []rPackage; 
}
