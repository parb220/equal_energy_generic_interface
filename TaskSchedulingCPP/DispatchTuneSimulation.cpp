#include <cmath>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <cstdio>
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

#include <time.h>

using namespace std; 

double CheckConvergency (std::vector<CSampleIDWeight> &samples, CEquiEnergyModel &model, int stage, int previous_stage,  double convergency_previous, double &average_consistency, double &std_consistency, double &LB_ESS, int posterior_type, int nGroup_NSE); 

vector<string> glob(const string &pattern);
void GetWeightedVarianceMatrix(CEquiEnergyModel &model, int stage, const std::vector<CSampleIDWeight> &); 
void AggregateScaleMatrix(CEquiEnergyModel &model, int stage); 
bool ScaleMatrixFileExist(CEquiEnergyModel &model, int stage); 
bool RenameScaleMatrixFile(CEquiEnergyModel &model, int p_stage, int c_stage); 
bool ScaleMatrixeFit(CEquiEnergyModel &model, int stage, int nNode, int nInitial, double lower_bound, double upper_bound );

std::vector<CSampleIDWeight> LoadSamplesFromFile(const string &file_name)
{
	CSampleIDWeight x;
        ifstream input_file;
        input_file.open(file_name.c_str(), ios::binary|ios::in);
        if (!input_file)
                return std::vector<CSampleIDWeight> (0); 
	std::vector<CSampleIDWeight> points; 
	while(!input_file.eof())
	{
        	read(input_file, &(x));
		points.push_back(x); 
	}
        input_file.close();
	return points; 
}

void DispatchTuneSimulation(int nNode, int nInitial, CEquiEnergyModel &model,const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag, int nGroup_NSE)
{
	double *sPackage = new double[N_MESSAGE], *rPackage = new double[N_MESSAGE];  
	MPI_Status status; 


	// diagnostic statistics
	// logMDD[i][0], constant of integration using ellipticals as proposals
	// logMDD[i][1], constant of integration using draws from the previous stage and importance.
	vector<vector<double> > logMDD(model.parameter->number_energy_stage+1, vector<double>(2, 0.0)); 
	vector<double> consistency(model.parameter->number_energy_stage, 0.0);
	vector<double> average_consistency(model.parameter->number_energy_stage, 0.0); 
	vector<double> std_consistency(model.parameter->number_energy_stage, 0.0);
	vector<double> LB_ESS(model.parameter->number_energy_stage, 0.0); 

	vector<CSampleIDWeight> samples;
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
	  	// cout << "DispatchTuneSimulation() - top of loop: stage=" << stage << " " << ctime(&rawtime) << endl;
		/////////////////////////////////////////////////////////////////////////////////
		// Highest + 1 stage 
		if (stage == model.parameter->highest_stage)
		{
			model.storage->InitializeBin(stage+1, model.current_sample.GetSize_Data());
			if (stage == model.parameter->number_energy_stage-1)
			  {
		                // draw highest stage + 1 sample
				time(&rawtime);
				// cout << "DispatchTuneSimulation() - drawing from prior: stage=" << stage+1 << " temperature: " << model.parameter->t[stage+1] << " " << ctime(&rawtime) << endl;

				// samples = samples of highest+1 stage
				samples = HighestPlus1Stage(nNode, nInitial, model);   // Sample from likelihood-heated only at the very high temperature 

				time(&rawtime);
				// cout << "DispatchTuneSimulation() - done drawing from prior: stage=" << stage+1 << " " << ctime(&rawtime) << endl;

				// compute log MDD
				time(&rawtime);			 
				// cout << "DispatchTuneSimulation() - computing MDD: stage=" << stage+1 << " " << ctime(&rawtime) << endl;

				/* get draws of the highest+1 stage
				samples.clear(); 
				if (!model.storage->DrawAllSample(stage+1, samples, unstructured, data_size))
				{
				    cerr << "DispatchTuneSimulation: error occurred when loading all samples.\n";
				    abort();
				}*/

				logMDD[stage+1][0] = LogMDD(samples, model, model.parameter->t[stage+1], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
				logMDD[stage+1][1] = 0.0;

				// open MDD file for output and write log MDD for stage + 1
				mdd_file.open(mdd_filename.c_str());
				if (!mdd_file.is_open())
				{
				    cerr << "DispatchTuneSimulation: Error opening " << mdd_filename << endl; 
				    abort(); 	
				}
				mdd_file << stage+1 << "\t" << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << endl; 
				// cout << setprecision(20) << "logMDD at stage " << stage+1 << ": " << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << endl; 

				time(&rawtime);
				// cout << "DispatchTuneSimulation() - done computing MDD: stage=" << stage << " " << ctime(&rawtime) << endl;
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
				// cout << setprecision(20) << "logMDD at stage " << mdd_stage << ": " << logMDD[mdd_stage][0] << "\t" << logMDD[mdd_stage][1] << endl;
			      }

			    // get samples of the previous stage 
			    samples.clear(); 
			    if (!model.storage->DrawAllSample(stage+1, samples, unstructured, data_size))
			      {
                		cerr << "DispatchTuneSimulation: error occurred when loading all samples.\n";
                       		abort();
			      }
			  }			

			// // debugging
			// ofstream out;
			// out.open("tmp.tmp");
			// for (unsigned int i=0; i < samples.size(); i++)
			//   out << model.log_prior_function(samples[i]) << endl; //samples[i].weight << '\t' << samples[i].reserved << endl;
			// out.close();
			// abort();
		}

		////////////////////////////////////////////////////////////////////////////////
		// Starting points
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - getting initial points: stage=" << stage << " " << ctime(&rawtime) << endl;
		vector<CSampleIDWeight> start_points(nInitial); 
		model.storage->InitializeBin(stage, model.current_sample.GetSize_Data()); 
		// if (model.storage->empty(stage+1) || !model.Initialize_WeightedSampling(samples, nInitial, stage+1, start_points)) // samples = samples of the previous stage 
		if (samples.empty() || !model.Initialize_WeightedSampling(samples, nInitial, stage+1, start_points)) // samples = samples of the previous stage 
		{
			for (int i=0; i<nInitial; i++)
				start_points[i] = mode; 
		} 
		model.storage->ClearDepositDrawHistory(stage+1); 
		model.storage->ClearStatus(stage+1); 
		model.storage->RestoreForFetch(stage+1); 

		string start_point_file; 
		// for (int i=0; i<nInitial; i++)
		// {
		convert.str(string());
        	convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << stage; // << "." << i;
        	start_point_file = model.parameter->storage_dir + convert.str();
        	output_file.open(start_point_file.c_str(), ios::binary|ios::out);
        	if (!output_file)
		{
               		cerr << "Error in writing to " << start_point_file << endl; 
			abort(); 	
		}
        	else
		{
			for (int i=0; i<(int)(start_points.size()); i++)
               			write(output_file, &(start_points[i]));
		}
        	output_file.close();
		//}
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - done getting initial points: stage=" << stage << " " << ctime(&rawtime) << endl;

		if (!ScaleMatrixFileExist(model, stage+1) || (RenameScaleMatrixFile(model, stage+1, stage) && !ScaleMatrixeFit(model, stage, nNode, nInitial, 0.234*0.8, 0.234*1.25)) )
		/////////////////////////////////////////////////////////////////////////////////
		// Send out tuning jobs when the existing scale matrix does not fit anymore
		{
			// Only need to tune on each slave node noce, because the results will be aggregated
			time(&rawtime);
			// cout << "DispatchTuneSimulation() - getting weighted variance matrix: stage=" << stage << " " << ctime(&rawtime) << endl;
		
			// samples = samples of the previous stage
			GetWeightedVarianceMatrix(model, stage, samples); 
		
			sPackage[LEVEL_INDEX] = stage; 
			sPackage[thin_INDEX] = model.parameter->thin; 
			sPackage[THIN_INDEX] = model.parameter->THIN; 
			sPackage[PEE_INDEX] = model.parameter->pee/(model.parameter->THIN/model.parameter->thin); 

			time(&rawtime);
			// cout << "DispatchTuneSimulation() - dispatching tuning: stage=" << stage << " " << ctime(&rawtime) << endl;

			for (int i=1; i<nNode; i++)
			{
				sPackage[GROUP_INDEX] = dw_uniform_int(nInitial); 
				sPackage[GROUP_NUMBER_INDEX] = 1; 
				MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD); 
			}

			for (int i=1; i<nNode; i++)
				MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, TUNE_TAG_BEFORE_SIMULATION, MPI_COMM_WORLD, &status);

			AggregateScaleMatrix(model, stage);
			time(&rawtime);
			// cout << "DispatchTuneSimulation() - done tuning: stage=" << stage << " " << ctime(&rawtime);
		}
	
		// Check Convergency
		if (stage == model.parameter->number_energy_stage-1 )
		{
			consistency[stage] = CheckConvergency(samples, model, stage, stage+1, 0.0, average_consistency[stage], std_consistency[stage], LB_ESS[stage], LIKELIHOOD_HEATED, nGroup_NSE);
		}
		else
		{
			consistency[stage] = CheckConvergency(samples, model, stage, stage+1, consistency[stage+1], average_consistency[stage], std_consistency[stage], LB_ESS[stage], LIKELIHOOD_HEATED, nGroup_NSE); 
		}
                cout << "Convergency Measure at Stage " << stage << ": " << setprecision(20) << consistency[stage] << "\t" << average_consistency[stage] << "\t" << std_consistency[stage]<< "\t" << LB_ESS[stage] << endl;  
		
		// logMDD using importance sampling	
		// samples = samples of the previous stage
		logMDD[stage][1] = LogMDD_Importance(samples, model, stage+1, stage, logMDD[stage+1][1], LIKELIHOOD_HEATED); 
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// simualtion
		//cout << "Simulation at " << stage << " for " << model.parameter->simulation_length << endl; 
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - dispatching simulation (" << model.parameter->simulation_length << "): stage=" << stage << " " << " temperature: " << model.parameter->t[stage] << " " << ctime(&rawtime) << endl;
		// samples = samples of the current stage
		samples = DispatchSimulation(nNode, nInitial, model, model.parameter->simulation_length, stage, SIMULATION_TAG) ;
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - done simulating (" << model.parameter->simulation_length << "): stage=" << stage << " " << ctime(&rawtime) << endl;
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// logMDD using bridge's method
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - computing MDD: stage=" << stage << " " << ctime(&rawtime) << endl;

	     	logMDD[stage][0] = LogMDD(samples, model, model.parameter->t[stage], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
		
		mdd_file << stage << "\t" << logMDD[stage][0] << "\t" << logMDD[stage][1] << endl;
		mdd_file.flush();
		cout << setprecision(20) << "logMDD at stage " << stage << ": " << logMDD[stage][0] << "\t" << logMDD[stage][1] << endl; 
		cout.flush();

		// to save space, remove stage+1 samples
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - deleting files: stage=" << stage << " " << ctime(&rawtime) << endl;
                if (save_space_flag  && stage > 0 ) //&& stage+1 < model.parameter->number_energy_stage-1 )
		{
			model.storage->ClearSample(stage+1);  
			convert.str(string());
                        convert << model.parameter->run_id << "/" << model.parameter->run_id << "*." << stage+1 << "*";;
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
		model.storage->ClearBin(stage+1); 
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - bottom of loop: stage=" << stage << " " << ctime(&rawtime) << endl;
	}

	mdd_file.close();

	delete []sPackage; 
	delete []rPackage; 
}

void GetWeightedVarianceMatrix(CEquiEnergyModel &model, int stage, const std::vector<CSampleIDWeight> &samples)
{
	// Calculate B for AdaptiveAfterSimulation_WeightedSampling_OnePass
	/*samples
	vector<CSampleIDWeight> samples;
        if (!model.storage->DrawAllSample(stage+1, samples) || samples.empty() )
        {
        	cerr << "DispatchTuneSimulation() : Error in loading samples of previous stage.\n";
        	abort();
        }*/

	// weight
	vector<double> log_weight = model.Reweight(samples, stage, stage+1);
       	double log_weight_sum = log_weight[0];
       	for (int i=1; i<(int)log_weight.size(); i++)
		log_weight_sum = AddLogs(log_weight_sum, log_weight[i]);
       	vector<double> weight(log_weight.size(), 0.0);
       	for (int i=0; i<(int)log_weight.size(); i++)
                weight[i] = exp(log_weight[i] - log_weight_sum);
			
	// block_scheme
	stringstream convert; 
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

void AggregateScaleMatrix(CEquiEnergyModel &model, int stage)
// Aggregate the scales obtained from computing nodes
{
	stringstream convert; 
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

	convert.str(string());
        convert << model.parameter->run_id << "/" << model.parameter->run_id << BMATRIX;
        string bmatrix_file_name = model.parameter->storage_dir  + convert.str();
	remove(bmatrix_file_name.c_str()); 
}

bool ScaleMatrixFileExist(CEquiEnergyModel &model, int stage)
{
	stringstream convert; 
        convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << stage;
        string block_file_name = model.parameter->storage_dir + convert.str();

	struct stat buffer;   
  	return (stat (block_file_name.c_str(), &buffer) == 0);			
}

bool RenameScaleMatrixFile(CEquiEnergyModel &model, int p_stage, int c_stage)
{
	stringstream convert; 
	convert.str(string()); 
        convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << p_stage;
        string p_name = model.parameter->storage_dir + convert.str();
	
	convert.str(string()); 
        convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << c_stage;
        string c_name = model.parameter->storage_dir + convert.str();

	return (rename(p_name.c_str(), c_name.c_str()) == 0); 	
}

bool ScaleMatrixeFit(CEquiEnergyModel &model, int stage, int nNode, int nInitial, double lower_bound, double upper_bound )
{
        double *sPackage = new double [N_MESSAGE];
        double *rPackage = new double [N_MESSAGE];
        sPackage[thin_INDEX] = 1;
        sPackage[THIN_INDEX] = 1;
        sPackage[LEVEL_INDEX] = stage;
        sPackage[PEE_INDEX] = 0;
       	sPackage[LENGTH_INDEX] = 0;
        sPackage[BURN_INDEX] = 1000;

        MPI_Status status;

	for (int i=1; i<nNode; i++)
	{
		sPackage[GROUP_INDEX] = dw_uniform_int(nInitial); 
		sPackage[GROUP_NUMBER_INDEX] = 1; 
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SCALE_MATRIX_FIT_TAG, MPI_COMM_WORLD); 
	}

	double avg_accpt_rate = 0.0; 
	for (int i=1; i<nNode; i++)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SCALE_MATRIX_FIT_TAG, MPI_COMM_WORLD, &status);
		avg_accpt_rate += rPackage[RETURN_INDEX_5]; 
	}
	avg_accpt_rate = avg_accpt_rate/(nNode-1); 

	delete [] sPackage; 
	delete [] rPackage; 
	cout << "Metropolis acceptance rate at stage " << stage << " using the scale matrix of the previous stage " << avg_accpt_rate << endl; 
	return (avg_accpt_rate >= lower_bound && avg_accpt_rate <= upper_bound); 
}


