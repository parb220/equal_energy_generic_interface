#include <cmath>
#include <algorithm>
#include <vector> 
#include <functional>
#include <ctime>
#include <fstream>
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "dw_dense_matrix.hpp"
#include "dw_math.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"

extern "C" {
	#include "dw_math.h"
	// #include "dw_switchio.h"
	#include "dw_rand.h"
}

using namespace std;

bool CEquiEnergyModel::MakeEquiEnergyJump(CSampleIDWeight &y_end, const CSampleIDWeight &y_initial)
{
	if(dw_uniform_rnd() <= parameter->pee && storage->DrawSample(energy_stage+1, storage->BinIndex(energy_stage+1,-y_initial.weight), y_end) ) // if a sample is successfully draw from bin
	{
		// calculate log_ratio in the current and the higher stages
		double log_ratio; 

		log_ratio = parameter->LogRatio_Stage(-(y_initial.weight), -(y_end.weight), energy_stage+1);
		log_ratio += parameter->LogRatio_Stage(-(y_end.weight), -(y_initial.weight), energy_stage); 
		if (log(dw_uniform_rnd()) <= log_ratio)
			return true; 
	}
	y_end = y_initial; 
	return false; 
}

void CEquiEnergyModel::SaveSampleToStorage(const CSampleIDWeight &sample)
{
        storage->DepositSample(energy_stage, storage->BinIndex(energy_stage, -sample.weight), sample);
}

void CEquiEnergyModel::Take_Sample_Just_Drawn_From_Storage(const CSampleIDWeight &x)
{
	current_sample = x; 
	current_sample.id = timer_when_started;
}

int CEquiEnergyModel::EE_Draw(int MH_thin)
{
	CSampleIDWeight x_new; 
	int new_sample_code = NO_JUMP; 
	double bounded_log_posterior_new; 

	if (MakeEquiEnergyJump(x_new, current_sample))
	{
		Take_Sample_Just_Drawn_From_Storage(x_new); 
		new_sample_code = EQUI_ENERGY_JUMP; 
	}
	else 
	{
		if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, MH_thin))
		{
			current_sample = x_new; 
			current_sample.id = timer_when_started;
			new_sample_code = METROPOLIS_JUMP; 
		}
	}
	
	return new_sample_code; 
}


double CEquiEnergyModel::BurnIn(int burn_in_length)
{
	CSampleIDWeight x_new; 
	int nJump =0; 
	double bounded_log_posterior_new; 
	for (int i=0; i<burn_in_length; i++)
	{
		if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, 1) )
		{
			current_sample = x_new;
			current_sample.id = timer_when_started;  
			nJump ++; 
		}
	}
	// cout << "MH Jump " << nJump << " out of " << burn_in_length << " in burning-in.\n"; 
	return (double)(nJump)/(double)(burn_in_length);   
}

void CEquiEnergyModel::Simulation_Prior(bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new(TDenseVector(current_sample.data.dim));   
	double bounded_log_posterior_new;
        bool if_write_file = false;
        ofstream output_file;
        if (!sample_file_name.empty() )
        {
                output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
                if (output_file)
                        if_write_file = true;
        }

	for (int i=0; i<parameter->simulation_length; i++)
	{
		do
		{
			DrawParametersFromPrior(x_new.data.vector); 
			x_new.DataChanged(); 
			log_posterior_function(x_new);
		} while (x_new.weight <= MINUS_INFINITY); 
		current_sample = x_new; 
		current_sample.id = timer_when_started;			
		if (if_storage)
               	      SaveSampleToStorage(current_sample);
                if (if_write_file)
                      write(output_file, &current_sample);
	}
	if (if_write_file)
                output_file.close();
}

void CEquiEnergyModel::Simulation_Within(bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new; 
	int nJump =0; 
	double bounded_log_posterior_new; 
	bool if_write_file = false; 
	ofstream output_file; 
	if (!sample_file_name.empty() )
	{
		output_file.open(sample_file_name.c_str(), ios::binary | ios::out); 
		if (output_file)
			if_write_file = true; 
	}

	for (int i=0; i<parameter->simulation_length; i++)
	{
		for (int j=0; j<parameter->thin; j++)
		{
			if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, parameter->THIN/parameter->thin) )
                	{
                        	current_sample = x_new;
                        	current_sample.id = timer_when_started;
                        	nJump ++;
                	}
		}
		
		if (if_storage)
			SaveSampleToStorage(current_sample);
		if (if_write_file)
			write(output_file, &current_sample); 
	}
	if (if_write_file)
		output_file.close(); 

	// cout << "MH Jump " << nJump << " out of " << parameter->simulation_length*parameter->thin << " in simulation.\n"; 
}

void CEquiEnergyModel::Simulation_Cross(bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new;

	int nEEJump=0, nMHJump=0; 
	bool if_write_file = false;
        ofstream output_file;
        if (!sample_file_name.empty() )
        {
                output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
                if (output_file)
                        if_write_file = true;
        }
	for (int i=0; i<parameter->simulation_length; i++)
	{
		for (int j=0; j<parameter->thin; j++)
		{
			int jump_code = EE_Draw(parameter->THIN/parameter->thin); 
			if (jump_code == EQUI_ENERGY_JUMP)
				nEEJump++; 
			else if (jump_code == METROPOLIS_JUMP)
				nMHJump++; 
		}	
	
		if (if_storage)
			SaveSampleToStorage(current_sample); 
		if (if_write_file)
			write(output_file, &current_sample); 
	}

	if (if_write_file)
		output_file.close(); 	

	// cout << "EE Jump " << nEEJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n"; 
	// cout << "MH Jump " << nMHJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n"; 
}

CEquiEnergyModel::CEquiEnergyModel() : 
gmm_mean(vector<TDenseVector>(0)), gmm_covariance_sqrt(vector<TDenseMatrix>(0)), 
gmm_covariance_sqrt_log_determinant(vector<double>(0)), gmm_covariance_sqrt_inverse(vector<TDenseMatrix>(0)), 
if_bounded(true), energy_stage(0), t_bound(1.0), current_sample(CSampleIDWeight()), timer_when_started(-1), metropolis(NULL), parameter(NULL), storage(NULL)
{}

CEquiEnergyModel::CEquiEnergyModel(bool _if_bounded, int eL, double _t, const CSampleIDWeight &_x, time_t _time, CMetropolis *_metropolis, CEESParameter *_parameter, CStorageHead *_storage) :
gmm_mean(vector<TDenseVector>(0)), gmm_covariance_sqrt(vector<TDenseMatrix>(0)), 
gmm_covariance_sqrt_log_determinant(vector<double>(0)), gmm_covariance_sqrt_inverse(vector<TDenseMatrix>(0)),
if_bounded(_if_bounded), energy_stage(eL), t_bound(_t), current_sample(_x), timer_when_started(_time), metropolis(_metropolis), parameter(_parameter), storage(_storage) 
{
}


