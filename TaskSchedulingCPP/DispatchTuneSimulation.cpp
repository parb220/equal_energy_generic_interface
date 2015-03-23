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

#include <time.h>

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

	stringstream convert; 
	ifstream input_file;
	ofstream output_file;

	// log MDD filename and stream
	string mdd_filename; 
	ofstream mdd_file;
	convert.str(string());
	convert << model.parameter->storage_dir << model.parameter->run_id << "/" << model.parameter->run_id << ".LogMDD.txt";
	mdd_filename=convert.str();

	time_t rawtime;

	for (int stage=model.parameter->highest_stage; stage>=model.parameter->lowest_stage; stage--)
	{
	  time(&rawtime);
	  cout << "DispatchTuneSimulation() - top of loop: stage=" << stage << " " << ctime(&rawtime) << endl;
		/////////////////////////////////////////////////////////////////////////////////

		/////////////////////////////////////////////////////////////////////////////////
                // this was the old code - hongwei suggested the code below
		// // Highest +1 stage 
		// if (stage == model.parameter->number_energy_stage-1)
		// {
		// 	HighestPlus1Stage(nNode, nInitial, model);
		// 	// HighestPlus1Stage_Prior(nNode, nInitial, model); // Sample from prior
		// 	/*posterior.clear(); 
		// 	if (!model.storage->DrawAllSample(stage+1, posterior, unstructured, data_size))
                //         {
                //                 cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
                //                 abort();
                //         }

		// 	// logMDD[stage+1][1]: logMDD using elliptical for draws from prior distribution
                //         logMDD[stage+1][0] = logMDD[stage+1][1] = logMDD[stage+1][2] = logMDD[stage+1][3] = LogMDD(posterior, model, model.parameter->t[stage+1], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
		// 	cout << setprecision(20) << "logMDD at stage " << stage+1 << ": " << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << "\t" << logMDD[stage+1][2] << "\t" << logMDD[stage+1][3] << endl; */
		// }

                ////////////////////////////////////////////////////////////////////////////////////
		// Highest + 1 stage 
		if (stage == model.parameter->highest_stage)
		{
			if (stage == model.parameter->number_energy_stage-1)
			  {
		                // draw highest stage + 1 sample
				time(&rawtime);
				cout << "DispatchTuneSimulation() - drawing from prior: stage=" << stage+1 << " temperature: " << model.parameter->t[stage+1] << " " << ctime(&rawtime) << endl;

				//HighestPlus1Stage_Prior(nNode, nInitial, model);   // Sample from prior
				HighestPlus1Stage(nNode, nInitial, model);    // Sample from prior with likelihood heated extremely

				time(&rawtime);
				cout << "DispatchTuneSimulation() - done drawing from prior: stage=" << stage << " " << ctime(&rawtime) << endl;

				// compute log MDD
				time(&rawtime);			 
				cout << "DispatchTuneSimulation() - computing MDD: stage=" << stage+1 << " " << ctime(&rawtime) << endl;

				// get posterior draws
				posterior.clear(); 
				if (!model.storage->DrawAllSample(stage+1, posterior, unstructured, data_size))
				  {
				    cerr << "DispatchTuneSimulation: error occurred when loading all samples.\n";
				    abort();
				  }

				//logMDD[stage+1][0] = LogMDD(posterior, model, model.parameter->t[stage+1], USE_TRUNCATED_POWER, PRIOR_ONLY);
				logMDD[stage+1][0] = LogMDD(posterior, model, model.parameter->t[stage+1], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
				logMDD[stage+1][1] = 0.0;

				// open MDD file for output and write log MDD for stage + 1
				mdd_file.open(mdd_filename.c_str());
				if (!mdd_file.is_open())
				  {
				    cerr << "DispatchTuneSimulation: Error opening " << mdd_filename << endl; 
				    abort(); 	
				  }
				mdd_file << stage+1 << "\t" << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << endl; 
				cout << setprecision(20) << "logMDD at stage " << stage+1 << ": " << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << endl; 

				time(&rawtime);
				cout << "DispatchTuneSimulation() - done computing MDD: stage=" << stage << " " << ctime(&rawtime) << endl;
			  }
			else
			  {
			    // read log MDD info
			    input_file.open(mdd_filename.c_str());
			    if (!input_file.is_open())
			      {
				cerr << "Error opening " << mdd_filename << endl; 
				abort(); 
			      }	
			    int mdd_stage=stage+2;
			    while (true)
			      if (!(input_file >> mdd_stage) || !(input_file  >> logMDD[mdd_stage][0] >> logMDD[mdd_stage][1])) break;
			    input_file.close();
			    if (mdd_stage > stage+1)
			      {
				cerr << "DispatchTuneSimulation: stage " << stage+1 << " not in " << mdd_filename << endl;
				abort();
			      }

			    // open MDD file for output and write log MDD
			    mdd_file.open(mdd_filename.c_str());
			    if (!mdd_file.is_open())
			      {
				cerr << "Error opening " << mdd_filename << endl; 
				abort(); 	
			      }
			    for (mdd_stage=logMDD.size()-1; mdd_stage >= stage+1; mdd_stage--)
			      {
				mdd_file << mdd_stage << "\t" << logMDD[mdd_stage][0] << "\t" << logMDD[mdd_stage][1] << endl;
				cout << setprecision(20) << "logMDD at stage " << mdd_stage << ": " << logMDD[mdd_stage][0] << "\t" << logMDD[mdd_stage][1] << endl;
			      }

			    // get posterior draws
			    posterior.clear(); 
			    if (!model.storage->DrawAllSample(stage+1, posterior, unstructured, data_size))
			      {
                		cerr << "DispatchTuneSimulation: error occurred when loading all samples.\n";
                       		abort();
			      }
			  }			

			// // debugging
			// ofstream out;
			// out.open("tmp.tmp");
			// for (unsigned int i=0; i < posterior.size(); i++)
			//   out << model.log_prior_function(posterior[i]) << endl; //posterior[i].weight << '\t' << posterior[i].reserved << endl;
			// out.close();
			// abort();

		}

		////////////////////////////////////////////////////////////////////////////////
		// Starting points
		time(&rawtime);
		cout << "DispatchTuneSimulation() - getting initial points: stage=" << stage << " " << ctime(&rawtime) << endl;
		vector<CSampleIDWeight> start_points(nInitial); 
		if (model.storage->empty(stage+1) || !model.Initialize_WeightedSampling(nInitial, stage+1, start_points))
		{
			for (int i=0; i<nInitial; i++)
				start_points[i] = mode; 
		} 
		model.storage->ClearDepositDrawHistory(stage+1); 
		model.storage->ClearStatus(stage+1); 
		model.storage->RestoreForFetch(stage+1); 

		string start_point_file; 
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
		time(&rawtime);
		cout << "DispatchTuneSimulation() - done getting initial points: stage=" << stage << " " << ctime(&rawtime) << endl;

		/////////////////////////////////////////////////////////////////////////////////
		// Send out tuning jobs
		// Only need to tune on each slave node noce, because the results will be aggregated
		if (stage != model.parameter->number_energy_stage)
		{
		  time(&rawtime);
		  cout << "DispatchTuneSimulation() - getting weighted variance matrix: stage=" << stage << " " << ctime(&rawtime) << endl;
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

		time(&rawtime);
		cout << "DispatchTuneSimulation() - dispatching tuning: stage=" << stage << " " << ctime(&rawtime) << endl;

		for (int i=1; i<nNode; i++)
		{
			sPackage[GROUP_INDEX] = dw_uniform_int(nInitial); 
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD); 
		}

		for (int i=1; i<nNode; i++)
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);


		time(&rawtime);
		cout << "DispatchTuneSimulation() - done tuning: stage=" << stage << " " << ctime(&rawtime);
		
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
		//cout << "Simulation at " << stage << " for " << model.parameter->simulation_length << endl; 
		time(&rawtime);
		cout << "DispatchTuneSimulation() - dispatching simulation (" << model.parameter->simulation_length << "): stage=" << stage << " " << " temperature: " << model.parameter->t[stage] << " " << ctime(&rawtime) << endl;
		DispatchSimulation(nNode, nInitial, model, model.parameter->simulation_length, stage, SIMULATION_TAG);
		time(&rawtime);
		cout << "DispatchTuneSimulation() - done simulating (" << model.parameter->simulation_length << "): stage=" << stage << " " << ctime(&rawtime) << endl;
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// logMDD
		// orginally, it was model.parameter->number_energy_stage - 10 because of problems with a certain matrix being singular - I'm
		// trying model.parameter->number_energy_stage -1 instead (both conditinals below

		time(&rawtime);
		cout << "DispatchTuneSimulation() - computing MDD: stage=" << stage << " " << ctime(&rawtime) << endl;

		// if (stage == model.parameter->number_energy_stage -10)
		// {
		// 	posterior.clear(); 
 		// 	if (!model.storage->DrawAllSample(stage+1, posterior, unstructured, data_size))
 		// 	{
 		// 		cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
 		// 		abort();
 		// 	}
 
		// 	// logMDD[stage+1][1]: logMDD using elliptical for draws from prior distribution
		// 	logMDD[stage+1][0] = logMDD[stage+1][1] = logMDD[stage+1][2] = logMDD[stage+1][3] = LogMDD(posterior, model, model.parameter->t[stage+1], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
		// 	cout << setprecision(20) << "logMDD at stage " << stage+1 << ": " << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << "\t" << logMDD[stage+1][2] << "\t" << logMDD[stage+1][3] << endl; 
		// }
		if (stage <= model.parameter->number_energy_stage - 1)
		{

		        // at model.parameter->highest_stage, posterior at stage+1 is loaded at top of loop, otherwise it
		        // is loaded in the previous loop below.
			logMDD[stage][1] = LogMDD_Importance(posterior, model, stage+1, stage, logMDD[stage+1][1], LIKELIHOOD_HEATED); 
			//proposal = posterior; 
			posterior.clear(); 
        		if (!model.storage->DrawAllSample(stage, posterior, unstructured, data_size))
        		{	
               			cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
               			abort();
        		}
	     		logMDD[stage][0] = LogMDD(posterior, model, model.parameter->t[stage], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
			//logMDD[stage][2] = LogMDD(proposal, posterior, model, stage+1, stage, logMDD[stage+1][0]); 
			//logMDD[stage][3] = LogMDD(proposal, posterior, model, stage+1, stage, logMDD[stage+1][3]); 
		
			mdd_file << stage << "\t" << logMDD[stage][0] << "\t" << logMDD[stage][1] << endl;
			cout << setprecision(20) << "logMDD at stage " << stage << ": " << logMDD[stage][0] << "\t" << logMDD[stage][1] << endl; 
		}

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
		time(&rawtime);
		cout << "DispatchTuneSimulation() - deleting files: stage=" << stage << " " << ctime(&rawtime) << endl;
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
		time(&rawtime);
		cout << "DispatchTuneSimulation() - bottom of loop: stage=" << stage << " " << ctime(&rawtime) << endl;
	
	}

	mdd_file.close();

	delete []sPackage; 
	delete []rPackage; 
}
