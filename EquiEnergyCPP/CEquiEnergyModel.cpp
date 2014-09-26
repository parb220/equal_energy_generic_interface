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

bool CEquiEnergyModel::Initialize_WeightedSampling(int N, int level_index, vector<CSampleIDWeight> &starters) const
{
        if (level_index < 1)
	        return false;

	if (starters.size() != N)
		starters.resize(N);

	// Get all previous level's samples
	vector<CSampleIDWeight> samples;  
	if (!storage->DrawAllSample(level_index, samples, true, current_sample.GetSize_Data()) || samples.size() == 0)
                return false;

	// Cumulative sum of importance weights
	vector<double> weight_sum(samples.size());
	double power=1.0/K(level_index-1)-1.0/K(level_index);
	
	// Accumulate logs
	weight_sum[0] = samples[0].weight*power;
	for (int i=1; i<(int)samples.size(); i++)
                weight_sum[i] = AddLogs(weight_sum[i-1], samples[i].weight*power);

	// Normalize and exponentiate
        double sum=weight_sum.back();
	for (int i=0; i<(int)samples.size(); i++)
		weight_sum[i] = exp(weight_sum[i] - sum); 

	for (int i=0; i<N; i++)
        {
	      double random_number = dw_uniform_rnd();
	      int position = std::lower_bound(weight_sum.begin(), weight_sum.end(), random_number)-weight_sum.begin();
	      starters[i] = samples[position];
	}
	return true; 
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

CEquiEnergyModel::CEquiEnergyModel(int rank, int _nNodes, TTimeSeries *_target, CEESParameter *_parameter) :
  target(_target), storage((CStorageHead*)NULL), parameter(_parameter), nParameters(target->NumberParameters()), 
  nNodes(_nNodes), energy_level(parameter->number_energy_level), current_sample(TDenseVector(nParameters),0,-1.0E300,false), 
  timer_when_started(time(NULL)), K(parameter->number_energy_level+1,-1.0), LogKernelIntegral(parameter->number_energy_level+1,0.0), 
  ring_counts(parameter->number_rings)
{
  storage = new CStorageHead(current_sample.GetSize_Data(), rank, parameter->run_id, parameter->storage_marker, parameter->storage_dir, parameter->number_energy_level);

  parameter->p_select = (parameter->expected_block_size <= 0.0) ? 1.0 : parameter->expected_block_size/(double)nParameters;
  parameter->Gn = 1 + (parameter->desired_G - 1)/(nNodes - 1);
  parameter->simulation_length=parameter->N * parameter->Gn * (nNodes - 1);
  parameter->pee_divided_psave = parameter->pee/parameter->p_save;
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
  Set the following fields
   (1) K(level)
   (2) scale (set to 1.0)
   (3) OrthonormalDirections
   (4) SqrtDiagonal
   (5) LogKernelIntegral(level)
   (6) ESS of importance sample from previous level
   (7) iImportanceSample (set to 0)

  If level is equal to parameter->number_energy_level or parameter->number_energy_level - 1,
  then the following fields are also set
   (8)  InitialOrthonormalDirections
   (9)  InitialSqrtDirections
   (10) InitialCenter

  In all cases, the initialization parameters are written to disk.
*/
void CEquiEnergyModel::SetupFromPreviousLevel(int level)
{
  energy_level=level;

  if (level < parameter->number_energy_level)
    ReadInitializationFile(level+1);

  if (level == parameter->number_energy_level)
    {
      K(level)=parameter->max_energy;
      LogKernelIntegral(level)=1.0;
      scale=1.0;
      InitialOrthonormalDirections=Identity(nParameters);
      InitialSqrtDiagonal=Ones(nParameters);
      InitialCenter=Zeros(nParameters);
      OrthonormalDirections=InitialOrthonormalDirections;
      SqrtDiagonal=InitialSqrtDiagonal;

      WriteInitializationFile();
    }
  else
    {
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
	  ReadInitializationFile(level+1);

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
	  LogKernelIntegral(level)=LogKernelIntegral(level+1) + ComputeLogKernelIntegral(log_kernel_target,log_density_proposal);
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
	  LogKernelIntegral(level)=LogKernelIntegral(level+1) + ComputeLogKernelIntegral(log_kernel_target,log_density_proposal);
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

      // scale
      scale=1.0;

      // write file
      WriteInitializationFile();
    }

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
	  if (dw_uniform_rnd() <= parameter->pee_divided_psave)  
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
	  current_sample.id = parameter->Gn*node;
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
  output_file << "Level: " << energy_level << endl;
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
    output_file << LogKernelIntegral(i) << ", ";
  output_file << endl;
  output_file << "Temperature:" << endl;
  for (int i=parameter->number_energy_level; i >= energy_level; i--)
    output_file << K(i) << ", ";
  output_file << endl;
    
  output_file.close();
}


/////////////////////////////////////////////////////////////////////////////////
void CEquiEnergyModel::WriteParameters(void)
{
  fstream out;
  OpenFile(out,true,"parameter");

  // machine readable form
  out << setprecision(12) << scientific;
  out << parameter->storage_dir << endl;
  out << parameter->storage_marker << endl;
  out << parameter->run_id << endl;
  out << parameter->number_energy_level << endl;
  out << parameter->pee << endl;
  out << parameter->highest_level << endl;
  out << parameter->lowest_level << endl;
  out << parameter->simulation_length << endl;
  out << parameter->expected_block_size << endl;
  out << parameter->p_select << endl;
  out << parameter->p_save << endl;
  out << parameter->pee_divided_psave << endl;
  out << parameter->tiny << endl;
  out << parameter->nImportanceSamples << endl;
  out << parameter->max_energy << endl;
  out << parameter->min_ess << endl;
  out << parameter->number_rings << endl;
  out << parameter->nu << endl;
  out << parameter->desired_G << endl;
  out << parameter->N << endl;
  out << parameter->Gn << endl;

  // human readable form
  out << setprecision(9) << fixed;
  out << "storage directory: " << parameter->storage_dir << endl;
  out << "storage marker: " << parameter->storage_marker << endl;
  out << "run id: " << parameter->run_id << endl;
  out << "number energy levels: " << parameter->number_energy_level << endl;
  out << "probability of ee jump: " << parameter->pee << endl;
  out << "highest level: " << parameter->highest_level << endl;
  out << "lowest level: " << parameter->lowest_level << endl;
  out << "simulation length: " << parameter->simulation_length << endl;
  out << "expected block size: " << parameter->expected_block_size << endl;
  out << "probability of selecting a parameter: " << parameter->p_select << endl;
  out << "probability of saving draw: " << parameter->p_save << endl;
  out << "additional wait to make ee jump " << parameter->pee_divided_psave << endl;
  out << "multiple for std. dev. of directions not selected: " << parameter->tiny << endl;
  out << "number of importance samples to make: " << parameter->nImportanceSamples << endl;
  out << "maximum temperature: " << parameter->max_energy << endl;
  out << "minimum ESS: " << parameter->min_ess << endl;
  out << "number rings: " << parameter->number_rings << endl;
  out << "degrees of freedom: " << parameter->nu << endl;
  out << "desired number of groups: " << parameter->desired_G << endl;
  out << "draws per group: " << parameter->N << endl;
  out << "groups per compute core: " << parameter->Gn << endl;

  out.close();
}

string CEquiEnergyModel::MakeFilename(const string &id, int level, int node)
{
  stringstream filename;
  filename << parameter->storage_dir << parameter->run_id << "/"  << parameter->run_id << '.' << id;
  if (level >= 0)
    {
      filename << '.' << level;
      if (node >= 0)
	filename << "." << node;
    }
  return filename.str();
}

void CEquiEnergyModel::OpenFile(fstream &file, bool output_file, const string &id, int level, int node)
{
  string filename=MakeFilename(id,level);
  file.open(filename.c_str(), output_file ? ios::out : ios::in );
  if (!file.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
}

/*
   The initialization file for the top level contains
     (1) level                                                     - energy_level
     (2) temperature                                               - K(model.energy_level)
     (3) log of the integral of the kernel                         - LogKernelIntegral(level)
     (4) scale                                                     - scale
     (5) othogonal transformation for independent univariate draws - InitialOrthonormalDirections
     (6) individual scale for independent univariate draws         - InitialSqrtDiagonal
     (7) center for initial distribuion                            - InitialCenter     

   The initialization file all other levels contains
     (1) level                                                     - energy_level
     (2) temperature                                               - K(model.energy_level)
     (3) log of the integral of the kernel                         - LogKernelIntegral(level)
     (4) scale                                                     - scale
     (5) othogonal transformation for independent univariate draws - OrthonormalDirections
     (6) individual scale for independent univariate draws         - SqrtDiagonal
     (7) center for initial distribuion                            - Center     
*/
void CEquiEnergyModel::WriteInitializationFile(void)
{
  fstream out;
  OpenFile(out,true,"initialization",energy_level);
  if (energy_level == parameter->number_energy_level)
    {
      out << scientific << setprecision(12);
      out << energy_level << endl;
      out << K(energy_level) << endl;
      out << LogKernelIntegral(energy_level) << endl;
      out << scale << endl;
      out << InitialOrthonormalDirections << endl;
      out << InitialSqrtDiagonal << endl;
      out << InitialCenter;
    }
  else
    {
      out << scientific << setprecision(12);
      out << energy_level << endl;
      out << K(energy_level) << endl;
      out << LogKernelIntegral(energy_level) << endl;
      out << scale << endl;
      out << OrthonormalDirections << endl;
      out << SqrtDiagonal << endl;
    }
  out.close();
}

/*
   Reads initialization file for level+1, if it exists, and level
*/
void CEquiEnergyModel::ReadInitializationFile(int level)
{
  fstream in;
  OpenFile(in,false,"initialization",level);
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
      in >> LogKernelIntegral.vector[level];
      in >> scale;
      InitialOrthonormalDirections.Resize(nParameters,nParameters);
      in >> InitialOrthonormalDirections;
      InitialSqrtDiagonal.Resize(nParameters);
      in >> InitialSqrtDiagonal;
      InitialCenter.Resize(nParameters);
      in >> InitialCenter;

      OrthonormalDirections=InitialOrthonormalDirections;
      SqrtDiagonal=InitialSqrtDiagonal;
    }
  else 
    {
      ReadInitializationFile(level+1);
      in >> K.vector[level];
      in >> LogKernelIntegral.vector[level];
      in >> scale;
      OrthonormalDirections.Resize(nParameters,nParameters);
      in >> OrthonormalDirections;
      SqrtDiagonal.Resize(nParameters);
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
  x.Resize(nParameters);
  double log_density=0.0;
  for (int j=nParameters-1; j >= 0; j--)
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
  for (int j=nParameters-1; j >= 0; j--)	
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
double CEquiEnergyModel::AnalyzeInitialDraws(void)
{
  if (energy_level != parameter->number_energy_level)
    {
      cerr << "AnalyzeInitialDraws(): energy_level must equal parameter->number_energy_level" << endl;
      abort();
    }

  int simulation_length=parameter->simulation_length;

  // read draws
  string filename;
  ifstream in;
  int node=0, count=0;
  TDenseVector x(nParameters), m_sum(nParameters,0.0);
  TDenseMatrix v_sum(nParameters,nParameters,0.0);
  while (true)
    {
      filename=MakeFilename("draws",energy_level,node);
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

  K(energy_level)=parameter->max_energy;
  scale=1.0;
  LogKernelIntegral(energy_level)=1.0;

  // write initialization file
  WriteInitializationFile();

  // compute ESS
  TDenseVector log_kernel_target(simulation_length), log_density_proposal(simulation_length);
  for (int i=simulation_length-1; i >= 0; i--)
    {	
      log_density_proposal(i)=InitialDraw(x);
      log_kernel_target(i)=target->LogPosterior(x)/K(energy_level);
    }
  return EffectiveSampleSize(log_kernel_target,log_density_proposal);
}
