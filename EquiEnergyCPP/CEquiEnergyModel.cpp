#include <cmath>
#include <algorithm>
#include <vector> 
#include <functional>
#include <ctime>
#include <fstream>
#include <iomanip> 

#include "CEESParameter.h"
#include "CStorageTail.h"
#include "CEquiEnergyModel.h"
#include "CStorageTail.h"
#include "dw_dense_matrix.hpp"
#include "dw_math.h"
#include "dw_rand.h"

using namespace std;

CEquiEnergyModel::CEquiEnergyModel(int rank, int _nNodes, TTimeSeries *_target, CEESParameter *_parameter) :
  target(_target), storage(new TStorage(target->NumberParameters(), _parameter->storage_dir, _parameter->run_id)), 
  parameter(_parameter), nParameters(target->NumberParameters()), nNodes(_nNodes), energy_level(parameter->number_energy_level), 
  current_sample(), K(parameter->number_energy_level+1,-1.0), LogKernelIntegral(parameter->number_energy_level+1,0.0), 
  ring_counts(parameter->number_rings), ESS(parameter->number_energy_level+1,-1.0)
{
  if (parameter->max_energy < 0.0) parameter->max_energy=10*target->NumberVariables()*target->NumberObservations();

  parameter->p_select = ((parameter->expected_block_size <= 0.0) || (parameter->expected_block_size >= nParameters)) ? 1.0 
    : parameter->expected_block_size/(double)nParameters;

  parameter->Gn = 1 + (parameter->desired_G - 1)/(nNodes - 1);

  if (parameter->p_save < parameter->pee)
    {
      parameter->N=(int)ceil((double)parameter->N*parameter->pee/parameter->p_save);
      parameter->p_save=parameter->pee;
    }
  parameter->pee_divided_psave = parameter->pee/parameter->p_save;

  parameter->simulation_length=parameter->N * parameter->Gn * (nNodes - 1);

  parameter->min_ess*=parameter->simulation_length;

  if (parameter->geometric >= 1.0)
    {
      parameter->geometric=pow(parameter->max_energy,-1.0/(double)(parameter->number_energy_level-1));
    }
}

void CEquiEnergyModel::Setup(int level, double Kplus)
{
  energy_level = level;
  scale = 1.0;
  if (level == parameter->number_energy_level-1)
    {
      ReadInitialDistribution();
      K(level)=parameter->max_energy;
      storage->SetupForInput(level+1,parameter->number_rings, K(level));
    }
  else
    {
      K(level+1) = Kplus;
      if (parameter->geometric >= 0)
	{
	  K(level)=parameter->geometric * K(level+1);
	  storage->SetupForInput(level+1, parameter->number_rings, K(level));
	}
      else
	K(level)=storage->SetupForInput(level+1, parameter->number_rings, Kplus, parameter->min_ess);
    }
  LogKernelIntegral(level)=storage->ComputeLogIntegralKernel();
  storage->ComputeWeightedVarianceDecomposition(OrthonormalDirections,SqrtDiagonal);
  //storage->ComputeVarianceDecomposition(OrthonormalDirections,SqrtDiagonal);
}

double CEquiEnergyModel::Tune(bool use_EE, double target_scale, int period, int max_period, bool verbose)
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
    max_period=64*period;

  // Beginning adaptive burn-in 
  while (period <= max_period)
    {
      // simulate approximately period/p_saved draws 
      if (use_EE)
	SimulateEE(false,G,N,false,false);
      else
	SimulateMH(false,G*N,false,false);

      // jump ratio
      previous_ratio = (double)nMHJumps/(double)nMHProposed;

      // output progress
      if (verbose)
	{
	  cout << setprecision(4) << "Tune(" << target_scale << ", " << period << ", " << max_period 
	       << ") acceptance = " << 100*previous_ratio << "  scale = " << scale << endl;
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

      if ((scale < 1.0e-20) || (scale > 1.0e20))
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
void CEquiEnergyModel::SimulateEE(bool store, int G, int N, bool start_from_current_sample, bool verbose)
{
  // verbose rules
  int verbose_factor = (G*N/10 < 1000) ? 1000 : N*G/10;

  // diagnostics
  int current_ring, new_ring;
  ring_counts.resize(storage->NumberRings());
  ring_changes_up.resize(storage->NumberRings());
  ring_changes_down.resize(storage->NumberRings());
  nEEJumps.resize(storage->NumberRings());
  nEEProposed.resize(storage->NumberRings());
  for (int i=0; i < ring_counts.size(); i++)
    ring_counts[i]=ring_changes_up[i]=ring_changes_down[i]=nEEJumps[i]=nEEProposed[i]=0;
  nMHJumps=nMHProposed=nRestarts=nSaved=0;

  // initial value
  if (!start_from_current_sample)
    {
      storage->ImportanceWeightedDraw(current_sample);
      nRestarts++;
    }
  current_ring=storage->Ring(current_sample.LogPosterior());

  // simulation runs
  TDraw x_new;
  for (int idx=0, g=0, verbose_count=verbose_factor, count=N; g < G; idx++)
    {       
      // verbose output
      if (verbose && (verbose_count  == 0))
	{
	  cout << setprecision(4) << "DSMH(" << node << "): " << " MH - " << nMHJumps << "/" << nMHProposed
	       << " = " << 100*(double)nMHJumps/(double)nMHProposed
	       << " - New starts = " << nRestarts << " - Number saved = " << nSaved << endl;
	  verbose_count=verbose_factor;
	}

      if (count == 0)
	{
	  if (++g < G)
	    {
	      storage->ImportanceWeightedDraw(current_sample);
	      current_ring=storage->Ring(current_sample.LogPosterior());
	      nRestarts++;
	      count=N;
	    }
	  continue;
	}

      // write previous draw
      if (dw_uniform_rnd() <= parameter->p_save)
	{
	  nSaved++;
	  count--;
	  verbose_count--;
	  current_sample.SetInfo(node*G + g, idx);
	  if (store) storage->AddDraw(current_sample); 

	  // make equi-energy jump
	  if (dw_uniform_rnd() <= parameter->pee_divided_psave)  
	    {
	      nEEProposed[current_ring]++;
	      storage->EqualWeightedDraw(x_new, current_ring);
	      if (energy_level < parameter->number_energy_level-1)
		{
		  if (log(dw_uniform_rnd()) <= (x_new.LogPosterior() - current_sample.LogPosterior())*(1.0/K(energy_level) - 1.0/K(energy_level+1)))
		    {
		      nEEJumps[current_ring]++;
		      x_new.SetLogKernel(x_new.LogPosterior()/K(energy_level) - LogKernelIntegral(energy_level));
		      current_sample = x_new; 			  
		    }
		}
	      else
		{
		  if (log(dw_uniform_rnd()) <= (x_new.LogPosterior() - current_sample.LogPosterior())/K(energy_level) 
		      - (x_new.LogKernel() - LogInitialDensity(current_sample.Parameters())))
		    {
		      nEEJumps[current_ring]++;
		      x_new.SetLogKernel(x_new.LogPosterior()/K(energy_level) - LogKernelIntegral(energy_level));
		      current_sample = x_new; 
		    }
		}

	      ring_counts[current_ring]++;
	      continue;
	    }   
	}

      // make MH jump
      nMHProposed++;
      GetMetropolisProposal(x_new);
      if (log(dw_uniform_rnd()) <= (x_new.LogPosterior() - current_sample.LogPosterior())/K(energy_level))
	{
	  nMHJumps++;
	  current_sample = x_new; 
	  new_ring=storage->Ring(current_sample.LogPosterior());
	  if (new_ring > current_ring)
	    {
	      ring_changes_up[current_ring]++;
	      current_ring=new_ring;
	    }
	  else if (new_ring < current_ring)
	    {
	      ring_changes_down[current_ring]++;
	      current_ring=new_ring;
	    }
	}
      ring_counts[current_ring]++;	
    }

  if (verbose)
    {
      cout << setprecision(4) << "DSMH(" << node << "): " << " MH - " << nMHJumps << "/" << nMHProposed 
	   << " = " << 100*(double)nMHJumps/(double)nMHProposed
	   << " - New starts = " << nRestarts << " - Number saved = " << nSaved << endl;
    }	
}

void CEquiEnergyModel::SimulateMH(bool store, int number_to_save, bool start_from_current_sample, bool verbose)
{
  // verbose rules
  int verbose_factor = (number_to_save/10 < 1000) ? 1000 : number_to_save/10;

  // diagnostics
  nMHJumps=nMHProposed=nSaved=0; 

  // initial value
  if (!start_from_current_sample)
    {
      storage->ImportanceWeightedDraw(current_sample);
    }

  // simulation runs
  TDraw x_new;
  for (int idx=0, verbose_count=verbose_factor; nSaved < number_to_save; idx++)
    {       
      // verbose output
      if (verbose && (verbose_count  == 0))
	{
	  cout << setprecision(4) << "MH(" << node << "): " << nMHJumps << "/" << nMHProposed << " = " << 100*(double)nMHJumps/(double)nMHProposed
	       << " -- Number saved = " << nSaved << endl;
	  verbose_count=verbose_factor;
	}

      // write previous draw?
      if (dw_uniform_rnd() <= parameter->p_save)
	{
	  nSaved++;
	  verbose_count--;
	  if (store)
	    {
	      current_sample.SetInfo(parameter->Gn*node, idx);
	      storage->AddDraw(current_sample);;
	    }
	}

      // make MH jump
      nMHProposed++;
      GetMetropolisProposal(x_new);
      if (log(dw_uniform_rnd()) <= (x_new.LogPosterior() - current_sample.LogPosterior())/K(energy_level))
	{
	  nMHJumps++;
	  current_sample = x_new; 
	}	
    }

  if (verbose)
    {
      cout << setprecision(4) << "MH(" << node << "): " << nMHJumps << "/" << nMHProposed << " = " << 100*(double)nMHJumps/(double)nMHProposed
	   << " -- Number saved = " << nSaved << endl;
    }	
}

double CEquiEnergyModel::GetMetropolisProposal(TDraw &new_x)
{
  int n=OrthonormalDirections.cols;
  TDenseVector v(n);

  for (int i=n-1; i >= 0; i--)
    v.vector[i]=scale*SqrtDiagonal(i)*((dw_uniform_rnd() <= parameter->p_select) ? dw_gaussian_rnd() : parameter->tiny * dw_gaussian_rnd());  
  v=current_sample.Parameters() + OrthonormalDirections*v;   
  double log_posterior=target->LogPosterior(v), 
    log_kernel=(energy_level < parameter->number_energy_level) ? log_posterior/K(energy_level) - LogKernelIntegral(energy_level) : LogInitialDensity(v);
  new_x.Set(v,log_posterior,log_kernel);

  return 0.0;
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void CEquiEnergyModel::WriteSimulationDiagnostic(int level, int node)
{
  fstream output_file;
  OpenFile(output_file,true,"Diagnostic",level,node);

  if (ring_counts.size() != storage->NumberRings())
    {
      cerr << "Incorrect number of rings?!?" << endl;
      abort();
    }

  output_file << setprecision(9) << scientific;
  output_file << energy_level << endl;
  output_file << K(energy_level) << endl;
  output_file << LogKernelIntegral(level) << endl;
  output_file << scale << endl;
  for (int i=0; i < ring_counts.size(); i++) output_file << ring_changes_up[i] << ' ';
  output_file << endl;
  for (int i=0; i < ring_counts.size(); i++) output_file << ring_changes_down[i] << ' ';
  output_file << endl;
  for (int i=0; i < ring_counts.size(); i++) output_file << ring_counts[i] << ' ';
  output_file << endl;
  for (int i=0; i < ring_counts.size(); i++) output_file << nEEJumps[i] << ' ';
  output_file << endl;
  for (int i=0; i < ring_counts.size(); i++) output_file << nEEProposed[i] << ' ';
  output_file  << endl;
  output_file << nMHJumps << endl;
  output_file << nMHProposed << endl;
  output_file << nRestarts << endl;
  output_file << nSaved << endl;
    
  output_file.close();
}

void CEquiEnergyModel::WriteSimulationDiagnostic(int level)
{
  int i_in;
  double x_in;
  ring_counts.resize(storage->NumberRings());
  ring_changes_up.resize(storage->NumberRings());
  ring_changes_down.resize(storage->NumberRings());
  nEEJumps.resize(storage->NumberRings());
  nEEProposed.resize(storage->NumberRings());
  for (int i=0; i < ring_counts.size(); i++)
    ring_counts[i]=ring_changes_up[i]=ring_changes_down[i]=nEEJumps[i]=nEEProposed[i]=0;
  nMHJumps=nMHProposed=nRestarts=nSaved=0; 
  int scale_in=0, K_in=0, integral_in=0;
  string filename;
  fstream file;
  int node=0;
  while (true)
    {
      filename=MakeFilename("Diagnostic",level,node);
      file.open(filename.c_str(),ios::in);
      if (!file.is_open()) break;
      file >> i_in;
      file >> x_in; if (fabs(x_in - K(energy_level)) > 1.0e-6) K_in++;
      file >> x_in; if (fabs(x_in - LogKernelIntegral(energy_level)) > 1.0e-6) integral_in++;
      file >> x_in; if (fabs(x_in - scale) > 1.0e-6) scale_in++;
      for (int i=0; i < ring_counts.size(); i++) { file >> i_in; ring_changes_up[i]+=i_in; }
      for (int i=0; i < ring_counts.size(); i++) { file >> i_in; ring_changes_down[i]+=i_in; }
      for (int i=0; i < ring_counts.size(); i++) { file >> i_in; ring_counts[i]+=i_in; }
      for (int i=0; i < ring_counts.size(); i++) { file >> i_in; nEEJumps[i]+=i_in; }
      for (int i=0; i < ring_counts.size(); i++) { file >> i_in; nEEProposed[i]+=i_in; }
      file >> i_in; nMHJumps+=i_in;
      file >> i_in; nMHProposed+=i_in;
      file >> i_in; nRestarts+=i_in;
      file >> i_in; nSaved+=i_in;
      file.close();
      remove(filename.c_str());
      node++;
    }
  if (node == 0)
    {
      cerr << "No diagnostic files - unable to open " << filename << endl;
      abort();
    }

  OpenFile(file,true,"Diagnostic",level);

  file << setprecision(8);
  file << "Stage: " << energy_level << endl;
  file << "Temperature: " << K(energy_level); if (K_in) file << " (" << K_in << ")" << endl; else file << endl;
  file << "Log of the integral of the kernel: " << LogKernelIntegral(energy_level); if (integral_in) file << " (" << integral_in << ")" << endl; else file << endl;
  file << "Scale: " << scale; if (scale_in) file << " (" << scale_in << ")" << endl; else file << endl;
  file << "Effective sample size: " << ESS(level) << endl;

  file << "Levels (in logs):" << endl;
  for (int i=1; i < storage->NumberRings(); i++) file << storage->LowerRingBoundary(i) << " ";
  file << endl;

  file << "Weighted striation probabilities:" << endl;
  for (int i=0; i < storage->NumberRings(); i++) file << "  " << 100*storage->ComputeWeightedStriationProbability(i) << " ";
  file << endl;

  file << "Ring changes up:" << endl;
  for (int i=0; i < storage->NumberRings(); i++) file << "  " << ring_changes_up[i];
  file << endl;

  file << "Ring changes down: " << endl;
  for (int i=0; i < storage->NumberRings(); i++) file << "  " << ring_changes_down[i];
  file << endl;

  file << "Ring counts:" << endl;
  for (int i=0; i < storage->NumberRings(); i++) file << "  " << ring_counts[i];
  file << endl;

  file << "Ring counts (fractions):" << endl;
  for (int i=0; i < storage->NumberRings(); i++) file << "  " << 100*(double)(ring_counts[i])/(double)(nMHProposed);
  file << endl;

  file << "Striated" << endl; 
  for (int i=0; i < storage->NumberRings(); i++) file << nEEJumps[i] << "/" << nEEProposed[i] << " = " << (nEEProposed[i] ? 100*(double)nEEJumps[i]/(double)nEEProposed[i] : 0.0) << endl;
  int eej=0, eep=0;
  for (int i=0; i < storage->NumberRings(); i++)
    {
      eej+=nEEJumps[i];
      eep+=nEEProposed[i];
    }
  file << "Total Striated: " << eej << "/" << eep << " = " << (eep ? 100*(double)eej/(double)eep : 0.0) << endl;

  file << "Metropolis-Hastings: " << nMHJumps << "/" << nMHProposed << " = " << (nMHProposed ? 100*(double)nMHJumps/(double)nMHProposed : 0.0) << endl;
  file << "Importance draw starts:" << nRestarts << endl;
  file << "Number saved: " << nSaved << endl;

  file << "Log kernel integral:" << endl;
  for (int i=parameter->number_energy_level; i >= energy_level; i--) file << LogKernelIntegral(i) << ", ";
  file << endl;

  file << "Temperature:" << endl; 
  for (int i=parameter->number_energy_level; i >= energy_level; i--) file << K(i) << ", ";
  file << endl;

  file << "Temperature ratio: " << endl;
  for (int i=parameter->number_energy_level-2; i >= energy_level; i--) file << K(i)/K(i+1) << ", ";
  file << endl;

  file << "Effective sample size:" << endl; 
  for (int i=parameter->number_energy_level; i >= energy_level; i--) file << ESS(i) << ", ";
  file << endl;

  file << endl << "SqrtDiagonal:" << endl;
  file << SqrtDiagonal << endl;
  file << "OrthonormalDirections:" << endl;
  file << OrthonormalDirections << endl;
    
  file.close();
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
  out << parameter->geometric << endl;

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
  out << "geometric progression (if positive): " << parameter->geometric << endl;

  out.close();
}

/////////////////////////////////////////////////////////////////////////////////
/// Initial distribution
/////////////////////////////////////////////////////////////////////////////////
void CEquiEnergyModel::WriteInitialDistribution(void)
{
  fstream out;
  OpenFile(out,true,"InitialDistribuiton");
  out << setprecision(9) << scientific;
  out << InitialCenter << endl;
  out << InitialSqrtDiagonal << endl;
  out << InitialOrthonormalDirections << endl;
  out.close();
}

void CEquiEnergyModel::ReadInitialDistribution(void)
{
  fstream in;
  OpenFile(in,false,"InitialDistribuiton");
  InitialCenter.Resize(nParameters);
  in >> InitialCenter;
  InitialSqrtDiagonal.Resize(nParameters);
  in >> InitialSqrtDiagonal;
  InitialOrthonormalDirections.Resize(nParameters,nParameters);
  in >> InitialOrthonormalDirections;
  in.close();
}

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
}

/////////////////////////////////////////////////////////////////////////////////
/// File creation utilities
/////////////////////////////////////////////////////////////////////////////////
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
  string filename=MakeFilename(id,level,node);
  file.open(filename.c_str(), output_file ? ios::out : ios::in );
  if (!file.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
}

