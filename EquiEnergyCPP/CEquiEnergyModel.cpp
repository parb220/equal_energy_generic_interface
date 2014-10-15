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

static double EffectiveSampleSize(double lambda, TDenseVector log_posterior, TDenseVector log_density)
{
  double sum=log_posterior.vector[0]*lambda - log_density.vector[0];
  for (int i=1; i < log_posterior.dim; i++)
    sum=AddLogs(sum,log_posterior.vector[i]*lambda - log_density.vector[i]);

  double weight, sum2=0;
  for (int i=0; i < log_posterior.dim; i++)
    {
      weight=exp(log_posterior.vector[i]*lambda - log_density.vector[i] - sum);
      sum2+=weight*weight;
    }

  return 1/sum2;
}


static double SetLambda(double lambda_min, double min_ess, TDenseVector log_posterior, TDenseVector log_density)
{
  // lambda = 1.0?
  if (EffectiveSampleSize(1.0,log_posterior,log_density) >= min_ess) return 1.0;

  // bisection
  double lambda, lambda_max=1.0;
  for (int i=0; i < 40; i++)
    {
      lambda=0.5*(lambda_max + lambda_min);
      if (EffectiveSampleSize(lambda,log_posterior,log_density) >= min_ess) 
	lambda_min=lambda;
      else
	lambda_max=lambda;
    }

  return 0.5*(lambda_max + lambda_min);
}

static double ComputeLogNormalization(double lambda, const TDenseVector &log_posterior, const TDenseVector &log_density)
{
  double sum=log_posterior.vector[0]*lambda - log_density.vector[0];
  for (int i=1; i < log_posterior.dim; i++)
    sum=AddLogs(sum,log_posterior.vector[i]*lambda - log_density.vector[i]);
  return sum;
}

static TDenseVector ComputeImportanceWeights(double lambda, const TDenseVector &log_posterior, const TDenseVector &log_density, double log_normalization)
{
  TDenseVector importance_weights(log_posterior.dim);
  for (int i=importance_weights.dim-1; i >= 0; i--)
    importance_weights.vector[i]=exp(log_posterior.vector[i]*lambda - log_density.vector[i] - log_normalization);
  return importance_weights;
}

static double ComputeEffectiveSampleSize(const TDenseVector &importance_weights)
{
  double sum2 = 0.0;
  for (int i=importance_weights.dim-1; i >= 0; i--)
    sum2+=importance_weights.vector[i]*importance_weights.vector[i];
  return 1.0/sum2;
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
CEquiEnergyModel::CEquiEnergyModel(int _node, int _nNodes, TTimeSeries *_target, CEESParameter *_parameter) :
  target(_target), storage(new TStorage(target->NumberParameters(), _parameter->storage_dir, _parameter->run_id, _node)), 
  parameter(_parameter), nParameters(target->NumberParameters()), nNodes(_nNodes), node(_node), stage(-1), current_sample(), 
  lambda(0), LogIntegralKernel(0), ESS(0), scale(1.0)
{
  if (parameter->min_lambda <= 0.0) parameter->min_lambda=1.0/(10.0*target->NumberVariables()*target->NumberObservations());

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
}

void CEquiEnergyModel::ComputeUnweightedMeanVariance(int stage, TDenseVector &mean, TDenseMatrix &variance)
{
  fstream file;
  storage->OpenFile(file,false,"draws",stage);

  int rows, cols;
  file.read((char*)(&rows),sizeof(int));
  file.read((char*)(&cols),sizeof(int));

  mean.Zeros(nParameters);
  variance.Zeros(nParameters,nParameters);
  TDenseVector v;
  for (int i=0; i < rows; i++)
    {
      file.seekg((2*sizeof(int)+sizeof(bool))+i*cols*sizeof(double),ios::beg);
      v.UniqueMemory(nParameters);
      file.read((char*)(v.vector),nParameters*sizeof(double));
      mean+=v;
      variance+=OuterProduct(v,v);
    }
  mean=(1.0/(double)rows)*mean;
  variance=(1.0/(double)rows)*variance;

  file.close();
}

void CEquiEnergyModel::ComputeWeightedMeanVariance(int stage, const TDenseVector &importance_weights, TDenseVector &mean, TDenseMatrix &variance)
{
  fstream file;
  storage->OpenFile(file,false,"draws",stage);

  int rows, cols;
  file.read((char*)(&rows),sizeof(int));
  file.read((char*)(&cols),sizeof(int));

  mean.Zeros(nParameters);
  variance.Zeros(nParameters,nParameters);
  TDenseVector v, w;
  for (int i=0; i < rows; i++)
    {
      file.seekg((2*sizeof(int)+sizeof(bool))+i*cols*sizeof(double),ios::beg);
      v.UniqueMemory(nParameters);
      file.read((char*)(v.vector),nParameters*sizeof(double));
      w=importance_weights.vector[i]*v;
      mean+=w;
      variance+=OuterProduct(w,v);
    }

  file.close();
}

void CEquiEnergyModel::WriteMHInfo(void)
{
  fstream file;
  storage->OpenFile(file,true,"MH-info",stage);

  file.write((char*)(&stage),sizeof(int));
  file.write((char*)(&scale),sizeof(double));

  InitialCenter.WriteBinary(file);
  InitialSqrtDiagonal.WriteBinary(file);
  InitialOrthonormalDirections.WriteBinary(file);

  lambda.WriteBinary(file);
  LogIntegralKernel.WriteBinary(file);
  ESS.WriteBinary(file);
  SqrtDiagonal.WriteBinary(file);
  OrthonormalDirections.WriteBinary(file);

  file.close();
}

void CEquiEnergyModel::ReadMHInfo(int Stage)
{
  fstream file;
  OpenFile(file,false,"MH-info",Stage);

  file.read((char*)(&stage),sizeof(int));
  if (stage != Stage)
    {
      cerr << "CEquiEnergyModel::ReadMHInfo(): internal stage index does not match filename" << endl;
      abort();
    }
  file.read((char*)(&scale),sizeof(double));

  InitialCenter.ReadBinary(file);
  InitialSqrtDiagonal.ReadBinary(file);
  InitialOrthonormalDirections.ReadBinary(file);

  lambda.ReadBinary(file);
  LogIntegralKernel.ReadBinary(file);
  ESS.ReadBinary(file);
  SqrtDiagonal.ReadBinary(file);
  OrthonormalDirections.ReadBinary(file);

  file.close();
}

void CEquiEnergyModel::SetupMasterNode(const TDenseVector &mean, const TDenseMatrix &variance)
{
  // sets stage
  stage=0;

  // sets InitialCenter, InitialSqrtDiagonal and InitialOrthonormalDirections
  InitialCenter=mean;
  Eig(InitialSqrtDiagonal,InitialOrthonormalDirections,variance);
  for (int i=InitialSqrtDiagonal.dim-1; i >= 0; i--)
    InitialSqrtDiagonal.vector[i]=sqrt(InitialSqrtDiagonal.vector[i]);

  // sets lambda
  lambda.Resize(1);
  lambda(0)=-1.0;

  // set LogIntegralKernel
  LogIntegralKernel.Resize(1);
  LogIntegralKernel(0)=0.0;

  // sets ESS
  ESS.Resize(1);
  ESS(0)=-1.0;

  // sets SqrtDiagonal and OrthonormalDirectons
  SqrtDiagonal=InitialSqrtDiagonal;
  OrthonormalDirections=InitialOrthonormalDirections;

  // sets scale
  scale=1.0;

  WriteMHInfo();
}

void CEquiEnergyModel::SetupMasterNode(int Stage)
{
  if (Stage <= 0)
    {
      cerr << "CEquiEnergyModel::SetupMasterNode(int Stage) - Stage must be greater than zero" << endl;
      abort(); 
    }

  // sets InitialCenter, InitialSqrtDiagonal, and InitialOrthonormalDirections 
  // sets initial segments of lambda, LogIntegralKernel, and ESS
  ReadMHInfo(Stage-1);
  TDenseVector log_posterior, log_density;
  storage->ReadInfo(Stage-1,log_posterior,log_density);

  // sets stage
  stage=Stage;

  // sets lambda(stage)
  if (stage == 1)
    lambda=Cat(lambda,parameter->min_lambda);
  else
    lambda=Cat(lambda,SetLambda(lambda(stage-1),parameter->min_ess,log_posterior,log_density));

  double log_normalization=ComputeLogNormalization(lambda(stage),log_posterior,log_density);
  TDenseVector importance_weights=ComputeImportanceWeights(lambda(stage),log_posterior,log_density,log_normalization);

  // sets LogIntegralKernel(stage)
  LogIntegralKernel=Cat(LogIntegralKernel,log_normalization - log(importance_weights.dim));

  // sets ESS(stage)
  ESS=Cat(ESS,ComputeEffectiveSampleSize(importance_weights));

  // sets SqrtDiagonal and OrthonormalDirections
  TDenseVector mean;
  TDenseMatrix variance;
  ComputeWeightedMeanVariance(stage-1,importance_weights,mean,variance);
  double min=1.0E300, max=-1.0E300;
  Eig(SqrtDiagonal,OrthonormalDirections,variance);
  for (int i=SqrtDiagonal.dim-1; i >= 0; i--)
    {
      if (SqrtDiagonal.vector[i] < min)	min=SqrtDiagonal.vector[i];
      if (SqrtDiagonal.vector[i] > max)	max=SqrtDiagonal.vector[i];
      SqrtDiagonal.vector[i]=sqrt(SqrtDiagonal.vector[i]);
    }
  if (max*1.0E-10 > min)
    {
      ComputeUnweightedMeanVariance(stage-1,mean,variance);
      Eig(SqrtDiagonal,OrthonormalDirections,variance);
      for (int i=SqrtDiagonal.dim-1; i >= 0; i--)
	SqrtDiagonal.vector[i]=sqrt(SqrtDiagonal.vector[i]);
    }

  // sets scale
  scale=1.0;

  WriteMHInfo();

  // setup storage
  storage->Setup(stage,parameter->number_striations,lambda.vector[stage]);
}

void CEquiEnergyModel::SetupComputeNode(int Stage)
{
  ReadMHInfo(Stage);
  storage->Setup(stage,parameter->number_striations,lambda.vector[stage]);
}

/*
   if use_EE is true, the following fields must be set:
     stage
     SqrtDiagonal
     OrthonormalDirections
     scale
     lambda(stage)                   if stage > 0
     LogIntegralKernel(stage)        if stage > 0
     InitialCenter                   if stage == 0
     InitialSqrtDiagonal             if stage == 0
     InitialOrthonormalDirections    if stage == 0

   if use_EE is false, the following fields must be set:
 
*/
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
      previous_ratio = nMHJumps.Sum()/nMHProposed.Sum();

      // output progress
      if (verbose)
	{
	  cout << setprecision(7) << "Tune(" << target_scale << ", " << period << ", " << max_period 
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
   The following fields must be set
*/
void CEquiEnergyModel::SimulateEE(bool store, int G, int N, bool start_from_current_sample, bool verbose)
{
  // verbose rules
  int verbose_factor = (G*N/10 < 1000) ? 1000 : N*G/10;

  // diagnostics
  int current_striation, new_striation;
  ResetDiagnostics(storage->NumberStriations());

  // initial value
  if (!start_from_current_sample)
    {
      storage->ImportanceWeightedDraw(current_sample);
      nStarts[current_striation=storage->Striation(current_sample.LogPosterior())]+=1;
    }
  else
    current_striation=storage->Striation(current_sample.LogPosterior());


  // simulation runs
  TDraw x_new;
  for (int idx=0, g=0, verbose_count=verbose_factor, count=N; g < G; idx++)
    {       
      // verbose output
      if (verbose && (verbose_count  == 0))
      	{
      	  cout << setprecision(7) << "DSMH(" << node << "): " << " MH - " << nMHJumps.Sum() << "/" << nMHProposed.Sum()
      	       << " = " << 100*nMHJumps.Sum()/nMHProposed.Sum()
      	       << " - New starts = " << nStarts.Sum() << " - Number saved = " << nSaved.Sum() << endl;
      	  verbose_count=verbose_factor;
      	}

      if (count == 0)
	{
	  if (++g < G)
	    {
	      storage->ImportanceWeightedDraw(current_sample);
	      nStarts[current_striation=storage->Striation(current_sample.LogPosterior())]+=1;
	      count=N;
	    }
	  continue;
	}

      // write previous draw
      if (dw_uniform_rnd() <= parameter->p_save)
	{
	  nSaved[current_striation]+=1;
	  count--;
	  verbose_count--;
	  current_sample.SetInfo(node*G + g, idx);
	  if (store) storage->AddDraw(current_sample); 

	  // make equi-energy jump
	  if (dw_uniform_rnd() <= parameter->pee_divided_psave)  
	    {
	      nEEProposed[current_striation]+=1;
	      storage->EqualWeightedDraw(x_new, current_striation);
	      if (stage > 1)
		{
		  if (log(dw_uniform_rnd()) <= (x_new.LogPosterior() - current_sample.LogPosterior())*(lambda(stage) - lambda(stage-1)))
		    {
		      nEEJumps[current_striation]+=1;
		      x_new.SetLogKernel(x_new.LogPosterior()*lambda(stage) - LogIntegralKernel(stage));
		      current_sample = x_new; 			  
		    }
		}
	      else
		{
		  if (log(dw_uniform_rnd()) <= (x_new.LogPosterior() - current_sample.LogPosterior())*lambda(stage) 
		      - (x_new.LogKernel() - LogInitialDensity(current_sample.Parameters())))
		    {
		      nEEJumps[current_striation]+=1;
		      x_new.SetLogKernel(x_new.LogPosterior()*lambda(stage) - LogIntegralKernel(stage));
		      current_sample = x_new; 
		    }
		}

	      striation_count[current_striation]+=1;
	      continue;
	    }   
	}

      // make MH jump
      nMHProposed[current_striation]+=1;
      GetMetropolisProposal(x_new);
      if (log(dw_uniform_rnd()) <= (x_new.LogPosterior() - current_sample.LogPosterior())*lambda(stage))
	{
	  nMHJumps[current_striation]+=1;
	  current_sample = x_new; 
	  new_striation=storage->Striation(current_sample.LogPosterior());
	  if (new_striation > current_striation)
	    {
	      striation_up[current_striation]+=1;
	      current_striation=new_striation;
	    }
	  else if (new_striation < current_striation)
	    {
	      striation_down[current_striation]+=1;
	      current_striation=new_striation;
	    }
	}
      striation_count[current_striation]+=1;	
    }

  if (verbose)
    {
      cout << setprecision(7) << "DSMH(" << node << "): " << " MH - " << nMHJumps.Sum() << "/" << nMHProposed.Sum()
	   << " = " << 100*nMHJumps.Sum()/nMHProposed.Sum()
	   << " - New starts = " << nStarts.Sum() << " - Number saved = " << nSaved.Sum() << endl;
    }	
}

void CEquiEnergyModel::SimulateMH(bool store, int number_to_save, bool start_from_current_sample, bool verbose)
{
  // verbose rules
  int verbose_factor = (number_to_save/10 < 1000) ? 1000 : number_to_save/10;

  // diagnostics
  ResetDiagnostics(1);

  // initial value
  if (!start_from_current_sample)
    {
      storage->ImportanceWeightedDraw(current_sample);
    }

  // simulation runs
  TDraw x_new;
  for (int idx=0, verbose_count=verbose_factor; nSaved(0) < number_to_save; idx++)
    {       
      // verbose output
      if (verbose && (verbose_count  == 0))
	{
	  cout << setprecision(7) << "MH(" << node << "): " << nMHJumps(0) << "/" << nMHProposed(0) << " = " << 100*(double)nMHJumps(0)/(double)nMHProposed(0) << " -- Number saved = " << nSaved(0) << endl;
	  verbose_count=verbose_factor;
	}

      // write previous draw?
      if (dw_uniform_rnd() <= parameter->p_save)
	{
	  nSaved(0)+=1;
	  verbose_count--;
	  if (store)
	    {
	      current_sample.SetInfo(parameter->Gn*node, idx);
	      storage->AddDraw(current_sample);;
	    }
	}

      // make MH jump
      nMHProposed(0)+=1;
      GetMetropolisProposal(x_new);
      if (log(dw_uniform_rnd()) <= (x_new.LogPosterior() - current_sample.LogPosterior())*lambda(stage))
	{
	  nMHJumps(0)+=1;
	  current_sample = x_new; 
	}	
    }

  if (verbose)
    {
      cout << setprecision(7) << "MH(" << node << "): " << nMHJumps(0) << "/" << nMHProposed(0) << " = " << 100*(double)nMHJumps(0)/(double)nMHProposed(0) << " -- Number saved = " << nSaved(0) << endl;
    }	
}

// The following fields must be set:
//   stage
//   SqrtDiagonal
//   OrthonormalDirections
//   scale
//   lambda(stage)
//   LogIntegralKernel(stage)
//   current_sample
//   parameter->p_select
//   parameter->tiny  
//
// If stage == 0, in addition to the above, the following fields must be set:
//     InitialCenter
//     InitialSqrtDiagonal
//     InitialOrthonormalDirections
//     parameter->nu
double CEquiEnergyModel::GetMetropolisProposal(TDraw &new_x)
{
  int n=OrthonormalDirections.cols;
  TDenseVector v(n);

  for (int i=n-1; i >= 0; i--)
    v.vector[i]=scale*SqrtDiagonal(i)*((dw_uniform_rnd() <= parameter->p_select) ? dw_gaussian_rnd() : parameter->tiny * dw_gaussian_rnd());  
  v=current_sample.Parameters() + OrthonormalDirections*v;   
  double log_posterior=target->LogPosterior(v), log_kernel=(stage > 0) ? log_posterior*lambda(stage) - LogIntegralKernel(stage) : LogInitialDensity(v);
  new_x.Set(v,log_posterior,log_kernel);

  return 0.0;
}

// The following fields must be set:
//   InitialCenter
//   InitialSqrtDiagonal
//   InitialOrthonormalDirections
//   parameter->nu
double CEquiEnergyModel::LogInitialDensity(const TDenseVector &x)
{
  double log_density=0.0;
  TDenseVector y=TransposeMultiply(InitialOrthonormalDirections,x-InitialCenter);
  for (int j=nParameters-1; j >= 0; j--)	
    log_density+=log(dw_tdistribution_pdf(y(j)/InitialSqrtDiagonal(j),parameter->nu)/InitialSqrtDiagonal(j));		
  return log_density;
}

// The following fields must be set:
//   InitialCenter
//   InitialSqrtDiagonal
//   InitialOrthonormalDirections
//   parameter->nu
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
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void CEquiEnergyModel::ResetDiagnostics(int number_striations)
{
  striation_count.Zeros(number_striations);
  striation_up.Zeros(number_striations);
  striation_down.Zeros(number_striations);
  nEEJumps.Zeros(number_striations);
  nEEProposed.Zeros(number_striations);
  nMHJumps.Zeros(number_striations);
  nMHProposed.Zeros(number_striations);
  nSaved.Zeros(number_striations);
  nStarts.Zeros(number_striations);
}

void CEquiEnergyModel::WriteDiagnosticComputeNode(void)
{
  fstream file;
  storage->OpenFile(file,true,"Diagnostic",stage,node);
  nEEJumps.WriteBinary(file);
  nEEProposed.WriteBinary(file);
  nMHJumps.WriteBinary(file);
  nMHProposed.WriteBinary(file);
  nSaved.WriteBinary(file);
  nStarts.WriteBinary(file);
  striation_count.WriteBinary(file);
  striation_up.WriteBinary(file);
  striation_down.WriteBinary(file);
  file.close();
}

void CEquiEnergyModel::WriteDiagnosticMasterNode(int stage)
{
  ResetDiagnostics(parameter->number_striations);
 
  TDenseVector x;
  string filename;
  fstream file;
  int node=1;
  while (true)
    {
      filename=storage->MakeFullFileName("Diagnostic",stage,node);
      file.open(filename.c_str(),ios::in|ios::binary);
      if (!file.is_open()) break;

      x.ReadBinary(file); nEEJumps+=x;
      x.ReadBinary(file); nEEProposed+=x;
      x.ReadBinary(file); nMHJumps+=x;
      x.ReadBinary(file); nMHProposed+=x;
      x.ReadBinary(file); nSaved+=x;
      x.ReadBinary(file); nStarts+=x;
      x.ReadBinary(file); striation_count+=x;
      x.ReadBinary(file); striation_up+=x;
      x.ReadBinary(file); striation_down+=x;

      file.close();
      remove(filename.c_str());
      node++;
    }
  if (node == 1)
    {
      cerr << "No diagnostic files - unable to open " << filename << endl;
      abort();
    }

  OpenFile(file,true,"Diagnostic",stage);

  file << setprecision(8);
  file << "Stage: " << stage << endl;
  file << "Lambda: " << lambda(stage) << endl;
  file << "Log of the integral of the kernel: " << LogIntegralKernel(stage) << endl;
  file << "Scale: " << scale << endl;
  file << "Effective sample size: " << ESS(stage) << endl;
  file << "Metropolis-Hastings acceptance: " << nMHJumps.Sum() << "/" << nMHProposed.Sum() << " = " << (nMHProposed.Sum() ? 100*nMHJumps.Sum()/nMHProposed.Sum() : 0.0) << endl;
  file << "Striated acceptance: " << nEEJumps.Sum() << "/" << nEEProposed.Sum() << " = " << (nEEProposed.Sum() ? 100*nEEJumps.Sum()/nEEProposed.Sum() : 0.0) << endl;
  file << "Number saved: " << nSaved.Sum() << endl;
  file << endl;

  file << "Weighted striation probabilities:" << endl;
  for (int i=0; i < storage->NumberStriations(); i++) file << "  " << 100*storage->WeightedStriationProbability(i);
  file << endl;
  file << "Striation counts (fractions):" << endl;
  for (int i=0; i < parameter->number_striations; i++) file << "  " << 100*striation_count(i)/nMHProposed.Sum();
  file << endl;
  file << "Striation counts (raw):" << endl << "  " << striation_count;
  file << "Striation changes up:" << endl << "  " << striation_up;
  file << "Striation changes down: " << endl << "  " << striation_down;
  file << endl;

  file << "Striated acceptance:" << endl; 
  for (int i=0; i < parameter->number_striations; i++) file << "  " << nEEJumps(i) << "/" << nEEProposed(i) << " = " << ((nEEProposed(i) > 0) ? 100*nEEJumps(i)/nEEProposed(i) : 0.0) << endl;
  file << endl;

  file << "Metropolis-Hastings acceptance:" << endl; 
  for (int i=0; i < parameter->number_striations; i++) file << "  " << nMHJumps(i) << "/" << nMHProposed(i) << " = " << ((nMHProposed(i) > 0) ? 100*nMHJumps(i)/nMHProposed(i) : 0.0) << endl;
  file << endl;

  file << "Importance draw starts:" << endl << "  " << nStarts << endl;
  file << "Number saved: " << endl << "  " << nSaved << endl;

  file << "Log kernel integral:" << endl << LogIntegralKernel(0);
  for (int i=1; i <= stage; i++) file << ", " <<  LogIntegralKernel(i);
  file << endl << endl;

  file << "lambda:" << endl << lambda(0); 
  for (int i=1; i <= stage; i++) file << ", " << lambda(i);
  file << endl << endl;

  if (stage > 1)
    {
      file << "lambda ratio: " << endl << lambda(2)/lambda(1);
      for (int i=2; i < stage; i++) file << ", " << lambda(i+1)/lambda(i);
      file << endl << endl;
    }

  file << "Effective sample size:" << endl << ESS(0); 
  for (int i=1; i <= stage; i++) file << ", " << ESS(i);
  file << endl << endl;

  file << "Striation boundaries (in logs):" << endl;
  for (int i=1; i < storage->NumberStriations(); i++) file << "  " << storage->LowerStriationBoundary(i);
  file << endl << endl;

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
  //out << parameter->number_energy_level << endl;
  out << parameter->pee << endl;
  out << parameter->first_stage << endl;
  out << parameter->simulation_length << endl;
  out << parameter->expected_block_size << endl;
  out << parameter->p_select << endl;
  out << parameter->p_save << endl;
  out << parameter->pee_divided_psave << endl;
  out << parameter->tiny << endl;
  out << parameter->nImportanceSamples << endl;
  out << parameter->min_lambda << endl;
  out << parameter->min_ess << endl;
  out << parameter->number_striations << endl;
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
  //out << "number energy levels: " << parameter->number_energy_level << endl;
  out << "probability of ee jump: " << parameter->pee << endl;
  out << "first stage: " << parameter->first_stage << endl;
  out << "simulation length: " << parameter->simulation_length << endl;
  out << "expected block size: " << parameter->expected_block_size << endl;
  out << "probability of selecting a parameter: " << parameter->p_select << endl;
  out << "probability of saving draw: " << parameter->p_save << endl;
  out << "additional wait to make ee jump " << parameter->pee_divided_psave << endl;
  out << "multiple for std. dev. of directions not selected: " << parameter->tiny << endl;
  out << "number of importance samples to make: " << parameter->nImportanceSamples << endl;
  out << "minimum lambda: " << parameter->min_lambda << endl;
  out << "minimum ESS: " << parameter->min_ess << endl;
  out << "number rings: " << parameter->number_striations << endl;
  out << "degrees of freedom: " << parameter->nu << endl;
  out << "desired number of groups: " << parameter->desired_G << endl;
  out << "draws per group: " << parameter->N << endl;
  out << "groups per compute core: " << parameter->Gn << endl;
  out << "geometric progression (if positive): " << parameter->geometric << endl;

  out.close();
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

