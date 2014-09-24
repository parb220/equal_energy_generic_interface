#include <cmath>
#include <algorithm>
#include <vector> 
#include <functional>
#include <ctime>
#include <fstream>
#include <iomanip> 

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
	/*if (energy_level == parameter->number_energy_level-1 && dw_uniform_rnd() <= parameter->pee && storage->DrawSample(energy_level+1, storage->BinIndex(energy_level+1,-y_initial.weight), y_end)  )
	{
		double log_ratio = GMM_LogRatio(y_initial, y_end); 
		if (log_ratio > MINUS_INFINITY)
		{
			log_ratio += parameter->LogRatio_Level(-y_end.weight, -y_initial.weight, energy_level); 
			if (log(dw_uniform_rnd()) <= log_ratio)
				return true; 
		}
		return true; 
	}
	// draw x_new from bin of the higher level of the same energy; 
	else if (energy_level<parameter->number_energy_level-1 && */ 
	if(dw_uniform_rnd() <= parameter->pee && storage->DrawSample(energy_level+1, storage->BinIndex(energy_level+1,-y_initial.weight), y_end) ) // if a sample is successfully draw from bin
	{
		// calculate log_ratio in the current and the higher levels
		double log_ratio = parameter->LogRatio_Level(-y_end.weight, -y_initial.weight, energy_level); 
		log_ratio += parameter->LogRatio_Level(-y_initial.weight, -y_end.weight, energy_level+1); 
		if (log(dw_uniform_rnd()) <= log_ratio)
			return true; 
	}
	y_end = y_initial; 
	return false; 
}

void CEquiEnergyModel::SaveSampleToStorage(const CSampleIDWeight &sample)
{
        storage->DepositSample(energy_level, storage->BinIndex(energy_level, -sample.weight), sample);
}

void CEquiEnergyModel::Take_Sample_Just_Drawn_From_Storage(const CSampleIDWeight &x)
{
	current_sample = x; 
	current_sample.id = (int)(time(NULL)-timer_when_started);
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
			current_sample.id = (int)(time(NULL)-timer_when_started);
			new_sample_code = METROPOLIS_JUMP; 
		}
	}
	
	return new_sample_code; 
}


double CEquiEnergyModel::BurnIn(int burn_in_length)
{
	CSampleIDWeight x_new; 
	int nJump =0; 
	double max_posterior = current_sample.weight, bounded_log_posterior_new; 
	for (int i=0; i<burn_in_length; i++)
	{
		if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, 1) )
		{
			current_sample = x_new;
			current_sample.id = (int)(time(NULL)-timer_when_started);  
			max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior; 
			nJump ++; 
		}
	}
	cout << "MH Jump " << nJump << " out of " << burn_in_length << " in burning-in.\n"; 
	return max_posterior;  
}

void CEquiEnergyModel::Simulation_Within(bool if_storage, const string &sample_file_name)
{
	// CSampleIDWeight x_new; 
	// int nJump =0; 
	// double bounded_log_posterior_new; 
	// bool if_write_file = false; 
	// ofstream output_file; 
	// if (!sample_file_name.empty() )
	// {
	// 	output_file.open(sample_file_name.c_str(), ios::binary | ios::out); 
	// 	if (output_file)
	// 		if_write_file = true; 
	// }

	// for (int i=0; i<parameter->simulation_length; i++)
	// {
	// 	for (int j=0; j<parameter->thin; j++)
	// 	{
	// 		if (metropolis->BlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, parameter->THIN/parameter->thin) )
        //         	{
        //                 	current_sample = x_new;
        //                 	current_sample.id = (int)(time(NULL)-timer_when_started);
        //                 	// max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior;
        //                 	nJump ++;
        //         	}
	// 	}
		
	// 	if (if_storage)
	// 		SaveSampleToStorage(current_sample);
	// 	if (if_write_file)
	// 		write(output_file, &current_sample); 
	// }
	// if (if_write_file)
	// 	output_file.close(); 

	// cout << "MH Jump " << nJump << " out of " << parameter->simulation_length*parameter->thin << " in simulation.\n"; 
}

void CEquiEnergyModel::Simulation_Cross(bool if_storage, const string &sample_file_name)
{
	// CSampleIDWeight x_new;

	// int nEEJump=0, nMHJump=0; 
	// bool if_write_file = false;
        // ofstream output_file;
        // if (!sample_file_name.empty() )
        // {
        //         output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
        //         if (output_file)
        //                 if_write_file = true;
        // }
	// for (int i=0; i<parameter->simulation_length; i++)
	// {
	// 	for (int j=0; j<parameter->thin; j++)
	// 	{
	// 		int jump_code = EE_Draw(parameter->THIN/parameter->thin); 
	// 		if (jump_code == EQUI_ENERGY_JUMP)
	// 			nEEJump++; 
	// 		else if (jump_code == METROPOLIS_JUMP)
	// 			nMHJump++; 
	// 		// if (jump_code == EQUI_ENERGY_JUMP || jump_code == METROPOLIS_JUMP)
	// 		//	max_posterior = max_posterior > current_sample.weight ? max_posterior : current_sample.weight; 
	// 	}	
	
	// 	if (if_storage)
	// 		SaveSampleToStorage(current_sample); 
	// 	if (if_write_file)
	// 		write(output_file, &current_sample); 
	// }

	// if (if_write_file)
	// 	output_file.close(); 	

	// cout << "EE Jump " << nEEJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n"; 
	// cout << "MH Jump " << nMHJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n"; 
}

CEquiEnergyModel::CEquiEnergyModel() : 
  gmm_mean(vector<TDenseVector>(0)), gmm_covariance_sqrt(vector<TDenseMatrix>(0)), 
  gmm_covariance_sqrt_log_determinant(vector<double>(0)), gmm_covariance_sqrt_inverse(vector<TDenseMatrix>(0)), 
  if_bounded(true), energy_level(0), t_bound(1.0), current_sample(CSampleIDWeight()), timer_when_started(time(NULL)), metropolis(NULL), parameter(NULL), storage(NULL),
  //K(0), IndependentDirections(0), scale(0)
K(0)
{};

CEquiEnergyModel::CEquiEnergyModel(bool _if_bounded, int eL, double _t, const CSampleIDWeight &_x, time_t _time, CMetropolis *_metropolis, CEESParameter *_parameter, CStorageHead *_storage) :
  gmm_mean(vector<TDenseVector>(0)), gmm_covariance_sqrt(vector<TDenseMatrix>(0)), 
  gmm_covariance_sqrt_log_determinant(vector<double>(0)), gmm_covariance_sqrt_inverse(vector<TDenseMatrix>(0)),
  if_bounded(_if_bounded), energy_level(eL), t_bound(_t), current_sample(_x), timer_when_started(_time), metropolis(_metropolis), parameter(_parameter), storage(_storage),
 //K(_parameter->number_energy_level+1), IndependentDirections(_parameter->number_energy_level+1), scale(_parameter->number_energy_level+1)
K(_parameter->number_energy_level+1)
{};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sampling and Tuning - added by DW 9/10/2014

double EffectiveSampleSize(double K1, double K2, vector<CSampleIDWeight> &samples)
{
  double power=1.0/K1 - 1.0/K2, sum=samples[0].weight*power, sum2=0.0, weight;

  for (int i=samples.size()-1; i > 0; i--)
    sum=AddLogs(sum,samples[i].weight*power);

  for (int i=samples.size()-1; i >= 0; i--)
    {
      weight=exp(samples[i].weight*power - sum);
      sum2+=weight*weight;
    }

  return 1.0/sum2;
}

double EffectiveSampleSize(TDenseVector log_kernel_target, TDenseVector log_density_proposal)
{
  if (log_kernel_target.dim != log_density_proposal.dim)
    {
      cerr << "EffectiveSampleSize(): log_kernel_target and log_density_proposal must be same length" << endl;
      abort();
    }

  double sum=log_kernel_target(0) - log_density_proposal(0), sum2=0.0, weight;
  for (int i=log_kernel_target.dim-1; i > 0; i--)
    sum=AddLogs(sum,log_kernel_target(i) - log_density_proposal(i));

  for (int i=log_kernel_target.dim-1; i >= 0; i--)
    {
      weight=exp(log_kernel_target(i) - log_density_proposal(i) - sum);
      sum2+=weight*weight;
    }

  return 1.0/sum2;
}

double ComputeLogKernelIntegral(TDenseVector log_kernel_target, TDenseVector log_density_proposal)
{
  if (log_kernel_target.dim != log_density_proposal.dim)
    {
      cerr << "LogKernelIntegral(): log_kernel_target and log_density_proposal must be same length" << endl;
      abort();
    }

  double sum=log_kernel_target(0) - log_density_proposal(0);
  for (int i=log_kernel_target.dim-1; i > 0; i--)
    sum=AddLogs(sum,log_kernel_target(i) - log_density_proposal(i));

  return sum - log(log_kernel_target.dim);
}


/*
   (1) Sets energy_level to level
   (2) Sets K(level) and IndependentDirections[level]
         Attempts to read these from file and otherwise computes them.  If 
         force_recompute is true, then no attempt is made to read their values 
         from disk.
   (3) Sets parameter->t[level] to K(level)
*/
void CEquiEnergyModel::SetupFromPreviousLevel(int level)
{
  if (level != energy_level)
    {
      cerr << "SetupFromPreviousLevel(): energy_level not equal to level" << endl;
      abort();
    }

  if (level == parameter->number_energy_level)
    {
      cerr << "Cannot setup highest level from previous level" << endl;
      abort();
    }

  // Get all previous level samples
  vector<CSampleIDWeight> samples;  
  if (!storage->DrawAllSample(level+1, samples, true, current_sample.GetSize_Data()) ) // || (samples.size() < parameter->simulation_length))
    {
      cerr << "Error obtaining all samples from level " << level+1 << endl;
      if (samples.size() < parameter->simulation_length) cerr << "Not enough samples - " << samples.size() << " : " << parameter->simulation_length << endl;
      abort();
    }

  TDenseVector log_kernel_target(samples.size()), log_density_proposal(samples.size()), log_posterior(samples.size());
  for (int i=samples.size()-1; i >= 0; i--)
    log_posterior(i)=samples[i].weight;

  if (level == parameter->number_energy_level-1)
    {
      energy_level++;
      ReadInitializationFile(level+1);
      energy_level--;

      K(level)=parameter->max_energy;
      for (int i=samples.size()-1; i >= 0; i--)	
	log_density_proposal(i)=LogInitialDensity(samples[i].data);
      log_kernel_target=(1.0/K(level))*log_posterior;
      ESS=EffectiveSampleSize(log_kernel_target,log_density_proposal);
      LogKernelIntegral(level)=ComputeLogKernelIntegral(log_kernel_target,log_density_proposal);
    }
  else if (level == 0)
    {
      K(level)=1.0;
      log_density_proposal=(1.0/K(level+1))*log_posterior;
      log_kernel_target=(1.0/K(level))*log_posterior;
      ESS=EffectiveSampleSize(log_kernel_target,log_density_proposal);
      LogKernelIntegral(level)=ComputeLogKernelIntegral(log_kernel_target,log_density_proposal);
    }
  else
    {
      double ess, K_new, K_max=K(level+1), K_min=(K_max/2 <= 1.0) ? 1.0 : K_max/2;

      // bracket min_ess
      do
	{
	  log_density_proposal=(1.0/K(level+1))*log_posterior;
	  log_kernel_target=(1.0/K_min)*log_posterior;
	  ess=EffectiveSampleSize(log_kernel_target,log_density_proposal);

	  cout << "K = " << K_min << "  ESS = " << ess << endl;

	  if (ess >= parameter->min_ess)  K_min=(K_min/2 <= 1.0) ? 1.0 : K_min/2;
	} 
      while ((K_min > 1) && (ess >= parameter->min_ess));

      // bracketed?
      if (ess < parameter->min_ess)
	{
	  for (int i=0; i < 10; i++)
	    {
	      K_new=0.5*(K_max + K_min);
	      log_density_proposal=(1.0/K(level+1))*log_posterior;
	      log_kernel_target=(1.0/K_new)*log_posterior;
	      ess=EffectiveSampleSize(log_kernel_target,log_density_proposal);

	      if (ess >= parameter->min_ess) 
		K_max=K_new;
	      else
		K_min=K_new;
	    }
	  K_new=0.5*(K_max + K_min);
	}
      else
	K_new=1.0;

      K(level)=K_new;
      log_density_proposal=(1.0/K(level+1))*log_posterior;
      log_kernel_target=(1.0/K(level))*log_posterior;
      ESS=EffectiveSampleSize(log_kernel_target,log_density_proposal);
      LogKernelIntegral(level)=ComputeLogKernelIntegral(log_kernel_target,log_density_proposal);
    }

  TDenseVector importance_weights(samples.size());
  TDenseMatrix variance;
  if (level < parameter->number_energy_level-1)
    {
      // Importance weights
      double power=(1.0/K(level)-1.0/K(level+1)), sum=samples[0].weight*power;
      for (int i=samples.size()-1; i > 0; i--)
  	sum=AddLogs(sum,samples[i].weight*power);
      for (int i=samples.size()-1; i >= 0; i--)
  	importance_weights(i)=exp(samples[i].weight*power - sum);

      // Compute variance matrix
      variance=OuterProduct(importance_weights(0)*samples[0].data,samples[0].data);
      for (int i=samples.size()-1; i > 0; i--)
	variance+=OuterProduct(importance_weights(i)*samples[i].data,samples[i].data);
      variance=0.5*(variance + Transpose(variance));
    }
  else
    {
      // TDenseVector log_initial_density(samples.size());
      // log_initial_density(0)=LogInitialDensity(samples[0].data);
      // double sum=samples[0].weight/K(level) - log_initial_density(0);
      // for (int i=samples.size()-1; i > 0; i--)
      // 	{
      // 	  log_initial_density(i)=LogInitialDensity(samples[i].data);
      // 	  sum=AddLogs(sum,samples[i].weight/K(level) - log_initial_density(i));
      // 	}
      // for (int i=samples.size()-1; i >= 0; i--)
      // 	importance_weights(i)=exp(samples[i].weight/K(level) - log_initial_density(i) - sum);

      // Compute variance matrix
      variance=OuterProduct(samples[0].data,samples[0].data);
      for (int i=samples.size()-1; i > 0; i--)
	variance+=OuterProduct(samples[i].data,samples[i].data);
      variance=0.5*(variance + Transpose(variance));
    }

  // Compute OrthonormalDirections and SqrtDiagonal
  Eig(SqrtDiagonal,OrthonormalDirections,variance);
  for (int j=SqrtDiagonal.dim-1; j >= 0; j--)
    SqrtDiagonal.vector[j]=sqrt(SqrtDiagonal.vector[j]);

  CreateInitializationFile(level,K(level),1.0,OrthonormalDirections,SqrtDiagonal);

  iImportanceSamples=0;
}

double CEquiEnergyModel::Tune(double target_scale, int period, int max_period, bool verbose)
{
  double previous_ratio; 
  double new_scale;
  double best_scale; 
  double low_scale; 
  double low_jump_ratio=-1.0; 
  double high_scale; 
  double high_jump_ratio=-1.0; 

  double log_mid=log(target_scale);
  double lower_bound=exp(log_mid/0.2);
  double upper_bound=exp(log_mid/5.0);

  int N=period;
  int G=1;

  if (max_period < period) 
    max_period=32*period;

  // Beginning adaptive burn-in 
  while (period <= max_period)
    {
      // simulate approximately period/p_saved draws 
      if (energy_level < parameter->number_energy_level)
	SimulateEE(false,string(),G,N,false,false);
      else
	SimulateMH(false,string(),G*N,false);

      // jump ratio
      previous_ratio = (double)nMHJumps/(double)nMHProposed;

      // output progress
      if (verbose)
	{
	  cout << setprecision(4) << "Tune(" << target_scale << ", " << period << ", " << max_period << ") acceptance = " << 100*previous_ratio << "  scale = " << scale << endl;  
	}

      // set new low or high bounds
      if (previous_ratio < target_scale)
	{
	  low_scale = scale; 
	  low_jump_ratio = previous_ratio; 
	}
      else 
	{
	  high_scale = scale; 
	  high_jump_ratio = previous_ratio; 
	}

      // new scale and best scale
      if (low_jump_ratio < 0.0)
	{
	  best_scale = scale; 
	  new_scale = (previous_ratio > upper_bound) ? 5.0*high_scale : (log_mid/log(previous_ratio))*high_scale;
	}
      else if (high_jump_ratio < 0.0)
	{
	  best_scale = scale; 
	  new_scale = (previous_ratio < lower_bound) ? 0.2*low_scale : (log_mid/log(previous_ratio))*low_scale;
	}
      else 
	{
	  best_scale = new_scale = ((target_scale-low_jump_ratio)*low_scale + (high_jump_ratio-target_scale)*high_scale)/(high_jump_ratio-low_jump_ratio);
	  period*=2; 
	  G*=2;
	  low_jump_ratio = -1.0; 
	  high_jump_ratio = -1.0; 
	}
			
      // reset scale
      scale = new_scale; 

      if ((scale < 1.0e-7) || (scale > 1.0e7))
	{
	  cout << "Unable to tune scale - returning" << endl;
	  return scale=best_scale;
	}
    }	  

  // set scale
  scale=best_scale;
  
  return scale;
}

/*
   Set reinitialize_factor to -1.0 to never reinitialize.
*/
void CEquiEnergyModel::SimulateEE(bool if_storage, const string &sample_file_name, int G, int N, bool start_from_current_sample, bool verbose)
{
  // setup file if needed
  bool if_write_file = false;
  ofstream output_file;
  if (!sample_file_name.empty() )
    {
      output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
      if (output_file)
	if_write_file = true;
    }

  // diagnostics
  int current_ring, new_ring;
  ring_counts.resize(storage->Number_Bin(energy_level+1));
  for (int i=0; i < ring_counts.size(); i++) ring_counts[i]=0;
  nEEJumps=nMHJumps=nEEProposed=nMHProposed=nRestarts=nSaved=ring_changes_up=ring_changes_down=0; 

  // initial value
  if (!start_from_current_sample)
    {
      ImportanceSamplePreviousLevel(current_sample);
      nRestarts++;
    }
  current_ring=storage->BinIndex(energy_level+1,-current_sample.weight);

  // simulation runs
  int verbose_factor = G*N/10;
  if (verbose_factor < 500) verbose_factor=500;
  CSampleIDWeight x_new;
  for (int g=0, verbose_count=verbose_factor, count=N; g < G; )
    {       
      // verbose output
      if (verbose && (verbose_count  == 0))
	{
	  cout << setprecision(4) << "EE: " << nEEJumps << "/" << nEEProposed << " = " << 100*(double)nEEJumps/(double)nEEProposed;
	  cout << " -- MH: " << nMHJumps << "/" << nMHProposed << " = " << 100*(double)nMHJumps/(double)nMHProposed;
	  cout << " -- Restarts = " << nRestarts << " -- Number saved = " << nSaved << endl;
	  verbose_count=verbose_factor;
	}

      if (count == 0)
	{
	  ImportanceSamplePreviousLevel(current_sample);
	  current_ring=storage->BinIndex(energy_level+1,-current_sample.weight);
	  g++;
	  count=N;
	  nRestarts++;
	}

      // write previous draw
      if (dw_uniform_rnd() <= parameter->p_save)
	{
	  nSaved++;
	  count--;
	  verbose_count--;
	  current_sample.id = node*G + g;
	  if (if_storage)      
	    SaveSampleToStorage(current_sample); 
	  if (if_write_file)
	    write(output_file, &current_sample);     

	  // make equi-energy jump
	  if (dw_uniform_rnd() <= parameter->pee)  
	    {
	      nEEProposed++;
	      if (storage->DrawSample(energy_level+1, current_ring, x_new))
		{
		  if (storage->BinIndex(energy_level+1,-current_sample.weight) != current_ring)
		    {
		      cerr << "DrawSample() did not draw from desired ring (ring=" << current_ring << ")" << endl;
		      abort();
		    }

		  if (energy_level < parameter->number_energy_level-1)
		    {
		      if (log(dw_uniform_rnd()) <= (x_new.weight - current_sample.weight)*(1.0/K(energy_level) - 1.0/K(energy_level+1)))
			{
			  nEEJumps++;
			  current_sample = x_new; 			  
			}
		    }
		  else
		    {
		      double log_initial_density_old=LogInitialDensity(current_sample.data), log_initial_density_new=LogInitialDensity(x_new.data);
		      if (log(dw_uniform_rnd()) <= (x_new.weight - current_sample.weight)/K(energy_level) - (log_initial_density_new - log_initial_density_old))
			{
			  nEEJumps++;
			  current_sample = x_new; 
			}
		    }
		}
	      else
		{
		  cerr << "Unable to obtain draw - level " << energy_level+1 << " - ring " << storage->BinIndex(energy_level+1,-current_sample.weight) << endl;
		  abort();
		}

	      ring_counts[current_ring]++;
	      continue;
	    }   
	}

      // make MH jump
      nMHProposed++;
      GetMetropolisProposal(x_new);
      if (log(dw_uniform_rnd()) <= (x_new.weight - current_sample.weight)/K(energy_level))
	{
	  nMHJumps++;
	  current_sample = x_new; 
	  new_ring=storage->BinIndex(energy_level+1,-current_sample.weight);
	  if (new_ring > current_ring)
	    {
	      ring_changes_up++;
	      current_ring=new_ring;
	    }
	  else if (new_ring < current_ring)
	    {
	      ring_changes_down++;
	      current_ring=new_ring;
	    }
	}
      ring_counts[current_ring]++;	
    }

  if (if_write_file)
    output_file.close(); 

  if (verbose)
    {
      cout << setprecision(4) << "EE: " << nEEJumps << "/" << nEEProposed << " = " << 100*(double)nEEJumps/(double)nEEProposed;
      cout << " -- MH: " << nMHJumps << "/" << nMHProposed << " = " << 100*(double)nMHJumps/(double)nMHProposed;
      cout << " -- Restarts = " << nRestarts << " -- Number saved = " << nSaved << endl;
    }	
}

void CEquiEnergyModel::SimulateMH(bool store_internal, const string &external_filename, int number_to_save, bool verbose)
{
  // setup file if needed
  bool store_external = false;
  ofstream output_file;
  if (!external_filename.empty() )
    {
      output_file.open(external_filename.c_str(), ios::out);
      if (output_file.is_open())
	{
	  output_file << scientific << setprecision(12);
	  store_external = true;
	}
    }

  // verbose rules
  int verbose_factor = (number_to_save/20 < 500) ? 500 : number_to_save/20;

  // diagnostics
  nMHJumps=nMHProposed=nSaved=0; 

  // simulation runs
  CSampleIDWeight x_new;
  for (int verbose_count=verbose_factor; nSaved < number_to_save; )
    {       
      // verbose output
      if (verbose && (verbose_count  == 0))
	{
	  cout << setprecision(4) << "MH: " << nMHJumps << "/" << nMHProposed << " = " << 100*(double)nMHJumps/(double)nMHProposed
	       << " -- Number saved = " << nSaved << endl;
	  verbose_count=verbose_factor;
	}

      // write previous draw?
      if (dw_uniform_rnd() <= parameter->p_save)
	{
	  nSaved++;
	  verbose_count--;
	  current_sample.id = parameter->G*node;
	  if (store_internal)
	    SaveSampleToStorage(current_sample); 
	  if (store_external)
	    output_file << current_sample.data;
	}

      // make MH jump
      nMHProposed++;
      GetMetropolisProposal(x_new);
      if (log(dw_uniform_rnd()) <= (x_new.weight - current_sample.weight)/K(energy_level))
	{
	  nMHJumps++;
	  current_sample = x_new; 
	}	
    }

  if (store_external)
    output_file.close(); 

  if (verbose)
    {
      cout << setprecision(4) << "MH: " << nMHJumps << "/" << nMHProposed << " = " << 100*(double)nMHJumps/(double)nMHProposed
	   << " -- Number saved = " << nSaved << endl;
    }	
}

double CEquiEnergyModel::GetMetropolisProposal(CSampleIDWeight &new_x)
{
  int n=OrthonormalDirections.cols;
  TDenseVector v(n);
  //double density=density_constant; 
  for (int i=n-1; i >= 0; i--)
    v.vector[i]=scale*SqrtDiagonal(i)*((dw_uniform_rnd() <= parameter->p_select) ? dw_gaussian_rnd() : parameter->tiny * dw_gaussian_rnd());     
  new_x.data=current_sample.data + OrthonormalDirections*v;
  new_x.weight=target->LogPosterior(new_x.data);
  new_x.calculated=true;

  return 0.0;
}

bool CEquiEnergyModel::ImportanceSamplePreviousLevel(CSampleIDWeight &x_new)
{
  if ((iImportanceSamples == 0) && !GetImportanceSamples()) return false;
  x_new=ImportanceSamples[--iImportanceSamples];
  return true;
}

bool CEquiEnergyModel::GetImportanceSamples(void)
{
  if (Initialize_WeightedSampling(parameter->nImportanceSamples,energy_level+1,ImportanceSamples))
    {
      iImportanceSamples=parameter->nImportanceSamples;
      return true;
    }
  else
    {
      cerr << "Unable to obtain importance samples from level " << energy_level+1 << endl;
      return false;
    }
}

void CEquiEnergyModel::WriteSimulationDiagnostic(int node)
{
  string filename=MakeFilename("Diagnostic",energy_level,node);
  ofstream output_file;
  output_file.open(filename.c_str(), ios::out);
  if (!output_file.is_open())
    {
      cerr << "Error opening " << filename << endl; 
      abort(); 	
    }

  output_file << setprecision(12) << scientific;
  output_file << K(energy_level) << endl;
  output_file << scale << endl;
  output_file << ring_changes_up << endl;
  output_file << ring_changes_down << endl;
  for (int i=0; i < ring_counts.size(); i++)
    output_file << ring_counts[i] << ' ';
  output_file << endl;
  output_file << nEEJumps << endl;
  output_file  << nEEProposed << endl;
  output_file << nMHJumps << endl;
  output_file << nMHProposed << endl;
  output_file << nRestarts << endl;
  output_file << nSaved << endl;
    
  output_file.close();
}

void CEquiEnergyModel::WriteSimulationDiagnostic(void)
{
  int i_in;
  double x_in;
  ring_counts.resize(parameter->number_rings);
  for (int i=0; i < ring_counts.size(); i++) ring_counts[i]=0;
  nEEJumps=nMHJumps=nEEProposed=nMHProposed=nRestarts=nSaved=ring_changes_up=ring_changes_down=0; 
  K(energy_level)=scale=0.0;
  string filename;
  ifstream in;
  int node=0;
  while (true)
    {
      filename=MakeFilename("Diagnostic",energy_level,node);
      in.open(filename.c_str(),ios::in);
      if (!in.is_open()) break;
      in >> x_in; K(energy_level)+=x_in;
      in >> x_in; scale+=x_in;
      in >> i_in; ring_changes_up+=i_in;
      in >> i_in; ring_changes_down+=i_in;
      for (int i=0; i < ring_counts.size(); i++)
	{
	  in >> i_in; ring_counts[i]+=i_in;
	}
      in >> i_in; nEEJumps+=i_in;
      in >> i_in; nEEProposed+=i_in;
      in >> i_in; nMHJumps+=i_in;
      in >> i_in; nMHProposed+=i_in;
      in >> i_in; nRestarts+=i_in;
      in >> i_in; nSaved+=i_in;
      in.close();
      node++;
    }
  if (node == 0)
    {
      cerr << "No diagnostic files - unable to open " << filename << endl;
      abort();
    }
  K(energy_level)/=node;
  scale/=node;

  ofstream output_file;
  filename=MakeFilename("Diagnostic",energy_level);
  output_file.open(filename.c_str(),ios::out);
  if (!output_file.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
  output_file << setprecision(8);
  output_file << "Temperature: " << K(energy_level) << endl;
  output_file << "Scale: " << scale << endl;
  output_file << "Effective sample size: " << ESS << endl;
  output_file << "Ring changes up: " << ring_changes_up << endl;
  output_file << "Ring changes down: " << ring_changes_down << endl;
  output_file << "Ring counts:" << endl;
  for (int i=0; i < ring_counts.size(); i++) 
    output_file << "  " << ring_counts[i];
  output_file << endl;
  output_file << "Ring counts:" << endl;
  for (int i=0; i < ring_counts.size(); i++) 
    output_file << "  " << 100*(double)(ring_counts[i])/(double)(nEEProposed + nMHProposed);
  output_file << endl;
  output_file << "Equi-energy: " << nEEJumps << "/" << nEEProposed << " = " << (nEEProposed ? 100*(double)nEEJumps/(double)nEEProposed : 0.0) << endl;
  output_file << "Metropolis-Hastings: " << nMHJumps << "/" << nMHProposed << " = " << (nMHProposed ? 100*(double)nMHJumps/(double)nMHProposed : 0.0) << endl;
  output_file << "Importance draw starts:" << nRestarts << endl;
  output_file << "Number saved: " << nSaved << endl;
  output_file << "Log kernel integral:" << endl;
  for (int i=parameter->number_energy_level; i >= energy_level; i--)
    output_file << LogKernelIntegral(i) << "  ";
  output_file << endl;
  output_file << "Temperature:" << endl;
  for (int i=parameter->number_energy_level; i >= energy_level; i--)
    output_file << K(i) << "  ";
  output_file << endl;
    
  output_file.close();
}


/////////////////////////////////////////////////////////////////////////////////
string CEquiEnergyModel::MakeFilename(const string &id, int level, int node)
{
  stringstream filename;
  filename << parameter->storage_dir << parameter->run_id << "/" << parameter->run_id << '.' << id << '.' << level << "." << node;
  return filename.str();
}

string CEquiEnergyModel::MakeFilename(const string &id, int level)
{
  stringstream filename;
  filename << parameter->storage_dir << parameter->run_id << "/" << parameter->run_id << '.' << id << '.' << level;
  return filename.str();
}

void CEquiEnergyModel::OpenFile(fstream &file, const string &id, int level, bool output_file)
{
  string filename=MakeFilename(id,level);
  file.open(filename, output_file ? ios::out : ios::in );
  if (!file.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
}

/*
   The initialization file contains
     (1) level                                                     - model.energy_level
     (2) temperature                                               - model.K(model.energy_level)
     (3) overall scale                                             - model.scale
     (4) othogonal transformation for independent univariate draws - model.OrthonormalDirections
     (5) individual scale for independent univariate draws         - model.SqrtDiagonal
*/
void CEquiEnergyModel::CreateInitializationFile(int level, double Te, double Sc, const TDenseMatrix &Or, const TDenseVector &Di)
{
  string filename=MakeFilename("initialization",level);
  ofstream out_file;
  out_file.open(filename.c_str(),ios::out);
  if (!out_file.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
  out_file << scientific << setprecision(12);
  out_file << level << endl;
  out_file << Te << endl;
  out_file << Sc << endl;
  out_file << Or << endl;
  out_file << Di;
  out_file.close();
}

/*
   The initialization file for the top level contains
     (1) level                                                     - model.energy_level
     (2) temperature                                               - model.K(model.energy_level)
     (3) othogonal transformation for independent univariate draws - model.InitialOrthonormalDirections
     (4) individual scale for independent univariate draws         - model.InitialSqrtDiagonal
     (5) center for initial distribuion                            - model.InitialCenter     
*/
void CEquiEnergyModel::CreateInitializationFile(int level, double Te, double Sc, const TDenseMatrix &Or, const TDenseVector &Di, const TDenseVector &Ce)
{
  if (level != parameter->number_energy_level)
    {
      cerr << "CreateInitializationFile(): This format initialization file only for top level" << endl;
      abort();
    }
  string filename=MakeFilename("initialization",level);
  ofstream out_file;
  out_file.open(filename.c_str(),ios::out);
  if (!out_file.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
  out_file << scientific << setprecision(12);
  out_file << level << endl;
  out_file << Te << endl;
  out_file << Sc << endl;
  out_file << Or << endl;
  out_file << Di << endl;
  out_file << Ce;
  out_file.close();
}

/*
   Reads initialization file and sets
     (1) model.energy_level - this is not yet done
     (2) model.K(model.energy_level)
     (3) model.scale
     (4) model.OrthonormalDirections
     (5) model.SqrtDiagonal
*/
void CEquiEnergyModel::ReadInitializationFile(int level)
{
  if (energy_level != level)
    {
      cerr << "energy_level not equal to level" << endl;
      abort();
    }
  string filename=MakeFilename("initialization",level);
  ifstream in;
  in.open(filename.c_str(),ios::in);
  if (!in.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
  int level_file;
  in >> level_file;
  if (level_file != level)
    {
      cerr << "level in file does not agree with filename " << endl;
      abort();
    }
  if (level == parameter->number_energy_level)
    {
      in >> K.vector[level];
      in >> scale;
      InitialOrthonormalDirections.Resize(parameter->nParameters,parameter->nParameters);
      in >> InitialOrthonormalDirections;
      InitialSqrtDiagonal.Resize(parameter->nParameters);
      in >> InitialSqrtDiagonal;
      InitialCenter.Resize(parameter->nParameters);
      in >> InitialCenter;

      SqrtDiagonal=InitialSqrtDiagonal;
      OrthonormalDirections=InitialOrthonormalDirections;
    }
  else 
    {
      if (level == parameter->number_energy_level-1)
	{
	  energy_level++;
	  ReadInitializationFile(level+1);
	  energy_level--;
	}

      in >> K.vector[level];
      in >> scale;
      OrthonormalDirections.Resize(parameter->nParameters,parameter->nParameters);
      in >> OrthonormalDirections;
      SqrtDiagonal.Resize(parameter->nParameters);
      in >> SqrtDiagonal;
    } 
  in.close();
  iImportanceSamples=0;
}

void CEquiEnergyModel::WriteScale(int level, int node, double s)
{
  string filename=MakeFilename("Scale",level,node);
  ofstream out;  
  out.open(filename.c_str(),ios::out);
  if (!out.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
  out << scientific << setprecision(12) << s;
  out.close();
}

void CEquiEnergyModel::WriteScale(int level, double s)
{
  string filename=MakeFilename("Scale",level);
  ofstream out;  
  out.open(filename.c_str(),ios::out);
  if (!out.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
  out << scientific << setprecision(12) << s;
  out.close();
}

double CEquiEnergyModel::ConsolidateScales(int level)
{
  double scale=0.0, scale_i;
  string filename;
  ifstream in;
  int node=0;
  while (true)
    {
      filename=MakeFilename("Scale",level,node);
      in.open(filename.c_str(),ios::in);
      if (!in.is_open()) break;
      in >> scale_i;
      in.close();
      scale+=scale_i;
      node++;
    }
  if (node == 0)
    {
      cerr << "No scale files - unable to open " << filename << endl;
      abort();
    }
  return scale/(double)node;
}

/////////////////////////////////////////////////////////////////////////////////
/// Initial distribution
/////////////////////////////////////////////////////////////////////////////////
double CEquiEnergyModel::InitialDraw(TDenseVector &x)
{
  x.Resize(parameter->nParameters);
  double log_density=0.0;
  for (int j=parameter->nParameters-1; j >= 0; j--)
    {
      x(j)=dw_tdistribution_rnd(parameter->nu);	
      log_density+=log(dw_tdistribution_pdf(x(j),parameter->nu)/InitialSqrtDiagonal(j));		
      x(j)*=InitialSqrtDiagonal(j);
    }
  x=InitialCenter + InitialOrthonormalDirections*x;
  return log_density;
}

double CEquiEnergyModel::LogInitialDensity(const TDenseVector &x)
{
  double log_density=0.0;
  TDenseVector y=TransposeMultiply(InitialOrthonormalDirections,x-InitialCenter);
  for (int j=parameter->nParameters-1; j >= 0; j--)	
    log_density+=log(dw_tdistribution_pdf(y(j)/InitialSqrtDiagonal(j),parameter->nu)/InitialSqrtDiagonal(j));		
  return log_density;

  //return target->LogPosterior(x)/K.vector[parameter->number_energy_level];
}


/*
   (1) Read draws from the top level
   (2) Compute InitialCenter, InitialOrthonormalDirections and InitialSqrtDiagonal
   (3) Returns effective sample size

   Note: K(level) and nu, the degress of freedom for the t-distribution, must be
   set upon entry.
*/
double CEquiEnergyModel::AnalyzeInitialDraws(int level)
{
  int nParameters=parameter->nParameters, simulation_length=parameter->simulation_length;

  // read draws
  string filename;
  ifstream in;
  int node=0, count=0;
  TDenseVector x(nParameters), m_sum(nParameters,0.0);
  TDenseMatrix v_sum(nParameters,nParameters,0.0);
  while (true)
    {
      filename=MakeFilename("draws",level,node);
      in.open(filename.c_str(),ios::in);
      if (!in.is_open()) break;

      while (true)
	{
	  try
	    {
	      in >> x;
	      m_sum+=x;
	      v_sum+=OuterProduct(x,x);
	      count++;
	    }
	  catch (dw_exception &e)
	    {
	      in.close();
	      break;
	    }
	}

      node++;
    }

  if (count == 0)
    {
      cerr << "Unable to read/process any draws file: " << filename << endl;
      abort();
    }

  // set InitialOrthonormalDirections, InitialSqrtDiagonal, and InitialCenter
  InitialCenter=(1.0/(double)count)*m_sum;
  Eig(InitialSqrtDiagonal,InitialOrthonormalDirections,(1.0/(double)count)*v_sum);
  for (int i=nParameters-1; i >= 0; i--)
    InitialSqrtDiagonal(i)=sqrt(InitialSqrtDiagonal(i));

  // write initialization file
  CreateInitializationFile(energy_level,parameter->max_energy,1.0,InitialOrthonormalDirections,InitialSqrtDiagonal,InitialCenter);

  // compute ESS
  TDenseVector log_kernel_target(simulation_length), log_density_proposal(simulation_length);
  for (int i=simulation_length-1; i >= 0; i--)
    {	
      log_density_proposal(i)=InitialDraw(x);
      log_kernel_target(i)=target->LogPosterior(x)/K(level);
    }
  return EffectiveSampleSize(log_kernel_target,log_density_proposal);
}
