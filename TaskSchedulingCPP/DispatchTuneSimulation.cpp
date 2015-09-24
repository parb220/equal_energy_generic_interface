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

using namespace std; 

vector<string> glob(const string &pattern);
void GetWeightedVarianceMatrix(CEquiEnergyModel &model, int stage, const std::vector<CSampleIDWeight> &); 
bool ReadScaleFromFile(const string &file_name, double &c); 
bool WriteScaleToFile(const string &file_name, double c); 
bool FileExist(const string &file_name); 
double ScaleFit(double*, double*, const int, CEquiEnergyModel &model, int stage, int nNode, int nInitial);
vector<int> StriationDistribution(const vector<CSampleIDWeight> &samples, const CEquiEnergyModel &model, int stage); 

void DispatchTuneSimulation(double *sPackage, double *rPackage, const int N_MESSAGE, int nNode, int nInitial, CEquiEnergyModel &model,const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag, int nGroup_NSE)
{
	// diagnostic statistics
	// logMDD[i][0], constant of integration using ellipticals as proposals
	// logMDD[i][1], constant of integration using draws from the previous stage and importance.
	vector<vector<double> > logMDD(model.parameter->number_energy_stage+1, vector<double>(2, 0.0)); 
	vector<double> consistency(model.parameter->number_energy_stage, 0.0), average_consistency(model.parameter->number_energy_stage, 0.0), std_consistency(model.parameter->number_energy_stage, 0.0), LB_ESS(model.parameter->number_energy_stage, 0.0); 
	vector<vector<int > > striation_distribution(model.parameter->number_energy_stage); 

	vector<CSampleIDWeight> samples;

	ifstream input_file;
	ofstream output_file, mdd_file;

	string mdd_filename = model.parameter->storage_dir + model.parameter->run_id + string("/") + model.parameter->run_id + string(".LogMDD.txt") ; 

	time_t rawtime;
	double alpha_0 = 0.234*0.8, alpha_1 = 0.234*1.25; 

	for (int stage=model.parameter->highest_stage; stage>=model.parameter->lowest_stage; stage--)
	{
	  	time(&rawtime);
	  	// cout << "DispatchTuneSimulation() - top of loop: stage=" << stage << " " << ctime(&rawtime) << endl;
		/////////////////////////////////////////////////////////////////////////////////
		// Highest + 1 stage 
		if (stage == model.parameter->highest_stage)
		{
			model.storage->InitializeBin(stage+1);
			if (stage == model.parameter->number_energy_stage-1)
			{
		                // draw highest stage + 1 sample
				time(&rawtime);
				// cout << "DispatchTuneSimulation() - drawing from prior: stage=" << stage+1 << " temperature: " << model.parameter->t[stage+1] << " " << ctime(&rawtime) << endl;

				// samples = samples of highest+1 stage
				samples = HighestPlus1Stage_Prior(sPackage, rPackage, N_MESSAGE, nNode, nInitial, model);   // Sample from prior

				time(&rawtime);
				// cout << "DispatchTuneSimulation() - done drawing from prior: stage=" << stage+1 << " " << ctime(&rawtime) << endl;

				// compute log MDD
				// cout << "DispatchTuneSimulation() - computing MDD: stage=" << stage+1 << " " << ctime(&rawtime) << endl;

				logMDD[stage+1][0] = LogMDD(samples, model, model.parameter->lambda[stage+1], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
				logMDD[stage+1][1] = 0.0;

				// open MDD file for output and write log MDD for stage + 1
				mdd_file.open(mdd_filename.c_str());
				if (!mdd_file.is_open())
				{
					cerr << "DispatchTuneSimulation: Error opening " << mdd_filename << endl; 
				    	abort(); 	
				}
				mdd_file << setprecision(20) << stage+1 << "\t" << logMDD[stage+1][0] << "\t" << logMDD[stage+1][1] << endl; 
				time(&rawtime);
				// cout << "DispatchTuneSimulation() - done computing MDD: stage=" << stage+1 << " " << ctime(&rawtime) << endl;
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
				{
			      		if (!(input_file >> mdd_stage) || !(input_file  >> logMDD[mdd_stage][0] >> logMDD[mdd_stage][1]))
						break;
				}
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
					mdd_file << setprecision(20) << mdd_stage << "\t" << logMDD[mdd_stage][0] << "\t" << logMDD[mdd_stage][1] << endl;

			    	// get samples of the previous stage 
			    	samples = model.storage->DrawAllSample(stage+1); 
			    	if (samples.empty())
			      	{
                			cerr << "DispatchTuneSimulation: error occurred when loading all samples.\n";
                       			abort();
			      	}
			  }			
		}
		// samples = samples of the previous stage
		// Check Convergency & logMDD with the importance weight method
		if (stage == model.parameter->number_energy_stage-1 )
			logMDD[stage][1] = consistency[stage] = CheckConvergency(samples, model, stage, stage+1, 0.0, average_consistency[stage], std_consistency[stage], LB_ESS[stage], LIKELIHOOD_HEATED, nGroup_NSE);
		else
			logMDD[stage][1] = consistency[stage] = CheckConvergency(samples, model, stage, stage+1, consistency[stage+1], average_consistency[stage], std_consistency[stage], LB_ESS[stage], LIKELIHOOD_HEATED, nGroup_NSE); 
                cout << "Convergency Measure at Stage " << stage << ": " << setprecision(20) << consistency[stage] << "\t" << average_consistency[stage] << "\t" << std_consistency[stage]<< "\t" << LB_ESS[stage] << endl;  
		
		////////////////////////////////////////////////////////////////////////////////
		// Starting points
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - getting initial points: stage=" << stage << " " << ctime(&rawtime) << endl;
		model.storage->InitializeBin(stage); 
		vector<CSampleIDWeight> start_points; 
		if (!samples.empty() )
			start_points = model.Initialize_WeightedSampling(samples, nInitial, stage+1);
		if (samples.empty() || start_points.empty())
		{
			start_points.resize(nInitial); 
			for (int i=0; i<nInitial; i++)
				start_points[i] = mode; 
		}

		string start_point_file = model.parameter->storage_dir + model.parameter->run_id + string("/") + model.parameter->run_id + START_POINT; 
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
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - done getting initial points: stage=" << stage << " " << ctime(&rawtime) << endl;


		/////////////////////// Tuning
		string block_file_name = model.parameter->storage_dir  + model.parameter->run_id + string("/") + model.parameter->run_id + BLOCK;
		if (!FileExist(block_file_name))
                        GetWeightedVarianceMatrix(model, stage, samples);
	
		double alpha = ScaleFit(sPackage, rPackage, N_MESSAGE, model, stage, nNode, nInitial); 
		int counter = 0; 
		while (alpha <= alpha_0 || alpha >= alpha_1) 
		{
			if (counter == 0)
				GetWeightedVarianceMatrix(model, stage, samples); 
			else 
			{
				if (log(alpha) <= 5.0*log(0.5*(alpha_0+alpha_1)))
					model.metropolis->SetScale(0.2); 
				else if (5.0*log(0.5*(alpha_0+alpha_1)) < log(alpha) && log(alpha) < 0.2*log(0.5*(alpha_0+alpha_1)))
					model.metropolis->SetScale(log(0.5*(alpha_0+alpha_1))/log(alpha)); 
				else if (log(alpha) >= 0.2*log(0.5*(alpha_0+alpha_1) ) )
					model.metropolis->SetScale(5.0); 

        			if (!model.metropolis->WriteBlocks(block_file_name))
        			{
                			cerr << "DispatchTuneSimulation() : Error in writing BMatrix file.\n";
                			abort();
        			}
			}
			
			counter ++; 

			alpha = ScaleFit(sPackage, rPackage, N_MESSAGE, model, stage, nNode, nInitial); 
		} 
		time(&rawtime);
	
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// simualtion
		// cout << "DispatchTuneSimulation() - dispatching simulation (" << model.parameter->simulation_length << "): stage=" << stage << " " << " temperature: " << model.parameter->t[stage] << " " << ctime(&rawtime) << endl;
		// samples = samples of the current stage
		samples = DispatchSimulation(sPackage, rPackage, N_MESSAGE, nNode, nInitial, model, model.parameter->simulation_length, stage, SIMULATION_TAG) ;
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - done simulating (" << model.parameter->simulation_length << "): stage=" << stage << " " << ctime(&rawtime) << endl;
	
		// Number of draws in each striation where the striation is defined at the pervious stage	
		striation_distribution[stage] = StriationDistribution(samples, model, stage+1);
		for (int i=0; i<(int)(striation_distribution[stage].size()); i++)
		{
			if (striation_distribution[stage][i])
				cout << "Number of sampels in striation " << i << " at stage " << stage << " " << striation_distribution[stage][i] << endl; 
		} 	
	
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// logMDD using bridge's method
		// cout << "DispatchTuneSimulation() - computing MDD: stage=" << stage << " " << ctime(&rawtime) << endl;
	     	logMDD[stage][0] = LogMDD(samples, model, model.parameter->lambda[stage], USE_TRUNCATED_POWER, LIKELIHOOD_HEATED);
		
		mdd_file << setprecision(20) << stage << "\t" << logMDD[stage][0] << "\t" << logMDD[stage][1] << endl;
		mdd_file.flush();
		cout << setprecision(20) << "logMDD at stage " << stage << ": " << logMDD[stage][0] << "\t" << logMDD[stage][1] << endl; 
		cout.flush();

		// to save space, remove stage+1 samples
		time(&rawtime);
		// cout << "DispatchTuneSimulation() - deleting files: stage=" << stage << " " << ctime(&rawtime) << endl;
                if (save_space_flag ) 
		{
			model.storage->ClearSample(stage+1);  
			stringstream convert; 
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
}

void GetWeightedVarianceMatrix(CEquiEnergyModel &model, int stage, const std::vector<CSampleIDWeight> &samples)
{
	// Calculate B for AdaptiveAfterSimulation_WeightedSampling_OnePass
	// weight
	vector<double> log_weight = model.Reweight(samples, stage, stage+1);
       	double log_weight_sum = log_weight[0];
       	for (int i=1; i<(int)log_weight.size(); i++)
		log_weight_sum = AddLogs(log_weight_sum, log_weight[i]);
       	vector<double> weight(log_weight.size(), 0.0);
       	for (int i=0; i<(int)log_weight.size(); i++)
                weight[i] = exp(log_weight[i] - log_weight_sum);
			
	// block_scheme
        string block_scheme_file_name = model.parameter->storage_dir  + model.parameter->run_id + string("/") + model.parameter->run_id + BLOCK_SCHEME;
	if (!model.metropolis->ReadBlockScheme(block_scheme_file_name)) 
		model.metropolis->SetBlockScheme(vector<TIndex>(1,TIndex(0,samples[0].data.dim-1))); 

	model.metropolis->GetBlockMatrix_WeightedSampling(samples, weight); 

	// Bmatrix file
	string bmatrix_file_name = model.parameter->storage_dir  + model.parameter->run_id  + string("/") + model.parameter->run_id + BLOCK; 
	if (!model.metropolis->WriteBlocks(bmatrix_file_name))
	{
		cerr << "DispatchTuneSimulation() : Error in writing BMatrix file.\n"; 
		abort(); 
	}
}

bool FileExist(const string &file_name)
{
	struct stat buffer;   
  	return (stat (file_name.c_str(), &buffer) == 0);			
}

bool ReadScaleFromFile(const string &file_name, double &c)
{
	ifstream iFile; 
	iFile.open(file_name.c_str(), iostream::in | iostream::binary); 
	if (!iFile)
		return false; 
	iFile.read((char *)&c, sizeof(double)); 
	iFile.close(); 
	return true; 
}

bool WriteScaleToFile(const string &file_name, double c)
{
	ofstream oFile; 
	oFile.open(file_name.c_str(), iostream::out | iostream::binary); 
	if (!oFile)
		return false; 
	oFile.write((char *)&c, sizeof(double)); 
	oFile.close(); 
	return true; 
}

double ScaleFit(double *sPackage, double *rPackage, const int N_MESSAGE, CEquiEnergyModel &model, int stage, int nNode, int nInitial)
{
        sPackage[THIN_INDEX] = 1;
        sPackage[LEVEL_INDEX] = stage;
        sPackage[PEE_INDEX] = model.parameter->pee;
       	sPackage[LENGTH_INDEX] = 500;
        sPackage[BURN_INDEX] = model.parameter->burn_in_length;

        MPI_Status status;

	for (int i=1; i<nNode; i++)
	{
		sPackage[GROUP_INDEX] = dw_uniform_int(nInitial); 
		sPackage[GROUP_NUMBER_INDEX] = 1; 
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, SCALE_MATRIX_FIT_TAG, MPI_COMM_WORLD); 
	}

	std::vector<int> nJump(2,0);
	for (int i=1; i<nNode; i++)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, SCALE_MATRIX_FIT_TAG, MPI_COMM_WORLD, &status);
		nJump[0] += rPackage[RETURN_INDEX_1]; // EE jump
                nJump[1] += rPackage[RETURN_INDEX_2]; // MH jump
	}
	double avg_accpt_rate = (double)nJump[1]/((nNode-1)*sPackage[LENGTH_INDEX]); 

	cout << "Metropolis acceptance rate at stage " << stage << " using the scale matrix of the previous stage " << avg_accpt_rate << endl; 
	return avg_accpt_rate; 
}

vector<int> StriationDistribution(const vector<CSampleIDWeight> &samples, const CEquiEnergyModel &model, int stage)
{
	int bin_index; 
	double log_posterior_stage; 
	vector<int> distribution(model.parameter->number_striation+1,0); 
	for (int i=0; i<(int)samples.size(); i++)
	{
		log_posterior_stage = samples[i].reserved*model.parameter->lambda[stage] + (samples[i].weight-samples[i].reserved); 
		bin_index = model.storage->BinIndex(stage,-log_posterior_stage); 
		distribution[bin_index]	++; 
	}
	return distribution; 
}

