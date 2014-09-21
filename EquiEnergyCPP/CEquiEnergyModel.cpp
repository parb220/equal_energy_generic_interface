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
                        	current_sample.id = (int)(time(NULL)-timer_when_started);
                        	// max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior;
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

	cout << "MH Jump " << nJump << " out of " << parameter->simulation_length*parameter->thin << " in simulation.\n"; 
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
			// if (jump_code == EQUI_ENERGY_JUMP || jump_code == METROPOLIS_JUMP)
			//	max_posterior = max_posterior > current_sample.weight ? max_posterior : current_sample.weight; 
		}	
	
		if (if_storage)
			SaveSampleToStorage(current_sample); 
		if (if_write_file)
			write(output_file, &current_sample); 
	}

	if (if_write_file)
		output_file.close(); 	

	cout << "EE Jump " << nEEJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n"; 
	cout << "MH Jump " << nMHJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n"; 
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

  if (sum2 < 0.0)
    {
      cerr << "Negative effective sample size: " << sum2 << endl;
      abort();
    }

  return 1.0/sum2;
}

double EffectiveSampleSize(double K1, TDenseVector log_posterior, TDenseVector log_initial_density)
{
  if (log_posterior.dim != log_initial_density.dim)
    {
      cerr << "EffectiveSampleSize(): log_posterior and log_initial_density must be same length" << endl;
      abort();
    }

  double sum=log_posterior(0)/K1 - log_initial_density(0), sum2=0.0, weight;
  for (int i=log_posterior.dim-1; i > 0; i--)
    sum=AddLogs(sum,log_posterior(i)/K1 - log_initial_density(i));

  for (int i=log_posterior.dim-1; i >= 0; i--)
    {
      weight=exp(log_posterior(i)/K1 - log_initial_density(i) - sum);
      sum2+=weight*weight;
    }

  return 1.0/sum2;
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
  if (level == parameter->number_energy_level)
    {
      cerr << "Cannot setup highest level from previous level" << endl;
      abort();
    }

  // Get all previous level samples
  cout << "SetupFromPreviousLevel(): calling DrawAllSample()\n";
  vector<CSampleIDWeight> samples;  
  if (!storage->DrawAllSample(level+1, samples, true, current_sample.GetSize_Data())) // || (samples.size() < parameter->simulation_length)) 
    {
      cerr << "Error obtaining all samples from level " << level+1 << endl;
      if (samples.size() < parameter->simulation_length)
	cerr << "Not enough samples - " << samples.size() << endl;
      abort();
    }

  if (level == parameter->number_energy_level-1)
    {
      K(level)=parameter->max_energy;
      K(level+1)=parameter->max_energy;
    }
  else if (level == 0)
    {
      K(level)=1.0;
    }
  else
    {
      // choose alpha so K(level+1) * alpha^(level+1) = 1
      double alpha=pow(K(level+1),-1.0/(level+1));

      K(level)=K(level+1)*alpha;
      double ess=EffectiveSampleSize(K(level),K(level+1),samples);

      if (ess > parameter->min_ess)
	{
	  double ess_new, K_new=K(level)*alpha;
	  int n=level;
	  while ((n-- > 0) && ((ess_new=EffectiveSampleSize(K_new,K(level+1),samples)) > parameter->min_ess))
	    {
	      K(level)=K_new;
	      ess=ess_new;
	      K_new*=alpha;
	    }
	}
    }
 
  // Importance weights
  TDenseVector importance_weights(samples.size());
  // if (level < parameter->number_energy_level-1)
  //   {
      double power=(1.0/K(level)-1.0/K(level+1)), sum=samples[0].weight*power;
      for (int i=samples.size()-1; i > 0; i--)
	sum=AddLogs(sum,samples[i].weight*power);
      for (int i=samples.size()-1; i >= 0; i--)
	importance_weights(i)=exp(samples[i].weight*power - sum);
  //   }
  // else
  //   {
  //     TDenseVector log_density(samples.size());
  //     log_density(0)=LogInitialDensity(samples[0].data);
  //     double sum=samples[0].weight/K(level) - log_density(0);
  //     for (int i=samples.size()-1; i > 0; i--)
  // 	{
  // 	  log_density(i)=LogInitialDensity(samples[i].data);
  // 	  sum=AddLogs(sum,samples[i].weight/K(level) - log_density(i));
  // 	}

  //     for (int i=samples.size()-1; i >= 0; i--)
  // 	importance_weights(i)=exp(samples[i].weight/K(level) - log_density(i) - sum);
  //   }

  // Compute variance matrix
  TDenseMatrix variance;
  variance=OuterProduct(importance_weights(0)*samples[0].data,samples[0].data);
  for (int i=samples.size()-1; i > 0; i--)
    variance+=OuterProduct(importance_weights(i)*samples[i].data,samples[i].data);
  variance=0.5*(variance + Transpose(variance));

  // Compute OrthonormalDirections and SqrtDiagonal
  Eig(SqrtDiagonal,OrthonormalDirections,variance);
  for (int j=SqrtDiagonal.dim-1; j >= 0; j--)
    SqrtDiagonal.vector[j]=sqrt(SqrtDiagonal.vector[j]);

  CreateInitializationFile(level,K(level),1.0,OrthonormalDirections,SqrtDiagonal);

  iImportanceSamples=0;

}

// /*
//    (1) Sets energy_level to level
//    (2) Sets K(level) and IndependentDirections[level]
//          Attempts to read these from file and otherwise computes them.  If 
//          force_recompute is true, then no attempt is made to read their values 
//          from disk.
//    (3) Sets parameter->t[level] to K(level)
// */
// bool CEquiEnergyModel::SetupLevel(int level, bool force_recompute)
// {
//   energy_level=level;

//   if (!force_recompute)
//     {
//       stringstream convert;
//       ifstream input_file;
//       convert << parameter->run_id << "/" << parameter->run_id << ".LevelSetup." << level;
//       string filename = parameter->storage_dir + convert.str();
//       input_file.open(filename.c_str(), ios::in);
//       if (input_file)
// 	{	  
// 	  try
// 	    {
// 	      input_file >> K(level);
// 	      input_file >> ess(level);
// 	      IndependentDirections[level].Resize(parameter->nParameters,parameter->nParameters);
// 	      input_file >> IndependentDirections[level];
// 	      input_file.close();
// 	      parameter->t[level]=K(level);
// 	      return true;
// 	    }
// 	  catch(dw_exception &e)
// 	    {
// 	      cerr << "Exception reading level setup file: " << filename << " : " << e.what() << endl;
// 	      input_file.close();
// 	      abort();
// 	    }
// 	}
//     }

//   // Get all previous level samples
//   vector<CSampleIDWeight> samples;  
//   if (!storage->DrawAllSample(level+1, samples, true, current_sample.GetSize_Data()) || (samples.size() < parameter->simulation_length)) 
//     {
//       cerr << "Error obtaining all samples from level " << level+1 << endl;
//       abort();
//     }

//   if (level == parameter->number_energy_level-1)
//     {
//       //K(level)=parameter->max_energy;
//       K(level+1)=parameter->t[level+1];
//       K(level)=parameter->t[level];
//       ess(level)=-1.0;
//     }
//   else if (level == 0)
//     {
//       K(level)=1.0;
//       ess(level)=EffectiveSampleSize(K(level),K(level+1),samples);
//     }
//   else
//     {
//       // choose alpha so K(level+1) * alpha^(level+1) = 1
//       double alpha=pow(K(level+1),-1.0/(level+1));

//       K(level)=K(level+1)*alpha;
//       ess(level)=EffectiveSampleSize(K(level),K(level+1),samples);

//       if (ess(level) < 0.0)
// 	{
// 	  cerr << "Negative effective sample size in SetupLevel(): " << ess(level) << endl;
// 	  abort();
// 	}

//       if (ess(level) > parameter->min_ess)
// 	{
// 	  double ess_new, K_new=K(level)*alpha;
// 	  int n=level;
// 	  while ((n-- > 0) && ((ess_new=EffectiveSampleSize(K_new,K(level+1),samples)) > parameter->min_ess))
// 	    {
// 	      K(level)=K_new;
// 	      ess(level)=ess_new;
// 	      K_new*=alpha;
// 	      if (ess(level) < 0.0)
// 		{
// 		  cerr << "Negative effective sample size in SetupLevel(): " << ess(level) << endl;
// 		  abort();
// 		}
// 	    }
// 	}
//     }

//   // Importance weights
//   TDenseVector importance_weights(samples.size());
//   double power=(1.0/K(level)-1.0/K(level+1));
//   double sum=samples[0].weight*power;
//   for (int i=samples.size()-1; i > 0; i--)
//     sum=AddLogs(sum,samples[i].weight*power);
//   for (int i=samples.size()-1; i >= 0; i--)
//     importance_weights(i)=exp(samples[i].weight*power - sum);

//   // Compute variance matrix
//   TDenseMatrix variance;
//   variance=OuterProduct(importance_weights(0)*samples[0].data,samples[0].data);
//   for (int i=samples.size()-1; i > 0; i--)
//     variance+=OuterProduct(importance_weights(i)*samples[i].data,samples[i].data);
//   variance=0.5*(variance + Transpose(variance));

//   // Compute independent directions
//   TDenseVector EigenValues;
//   Eig(EigenValues,IndependentDirections[level],variance);
//   for (int j=EigenValues.dim-1; j >= 0; j--)
//     {
//       double factor=sqrt(EigenValues(j));
//       for (int i=EigenValues.dim-1; i >= 0; i--)
// 	IndependentDirections[level](i,j)*=factor;
//     }

//   // Write K(level) and independent directions
//   stringstream convert;
//   ofstream output_file;
//   convert << parameter->run_id << "/" << parameter->run_id << ".LevelSetup." << level;
//   string filename = parameter->storage_dir + convert.str();
//   output_file.open(filename.c_str(), ios::out);
//   if (!output_file)
//     {
//       cerr << "Error in opening " << filename << endl; 
//       return false; 	
//     }
//   else
//     {
//       output_file << K(level) << endl;
//       output_file << ess(level) << endl;
//       output_file << IndependentDirections[level];
//     }

//   parameter->t[level]=K(level);

//   return true; 

// }


// /*
//   (0) SetupLevel() must be called prior to SetScale()
//   (1) Sets scale(energy_level)
//         (a) if new_scale > 0, scale(energy_level) = new_scale
//         (b) if new_scale <= 0, reads scale(energy_level) from file
//   (2) Computes density_constant for metropolis jumps (not currently implemented)
//   (3) Sets ImportanceSamples

//   After a successful call to SetScales()
// */ 
// bool CEquiEnergyModel::SetScale(double new_scale)
// {
//   // sets scale
//   if (new_scale <= 0.0)
//     {
//       stringstream convert;
//       convert << parameter->run_id << "/" << parameter->run_id << ".Scale." << energy_level;
//       string filename = parameter->storage_dir + convert.str();
//       ifstream input_file;
//       input_file.open(filename.c_str(), ios::in);
//       if (!input_file)
// 	{
// 	  cerr << "Error in opening " << filename << endl; 
// 	  abort();	
// 	}
//       else
// 	{
// 	  try
// 	    {
// 	      input_file >> scale(energy_level);
// 	    }
// 	  catch(dw_exception &e)
// 	    {
// 	      cerr << "Exception reading scale from " << filename << endl << e.what() << endl;
// 	      abort();
// 	    }
// 	}
//       input_file.close();
//     }
//   else
//     scale(energy_level)=new_scale;

//   // sets density_constant
//   density_constant=0.0;

//   // gets importance sample from previous level
//   return GetImportanceSamples();
// }

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

  int reinitialize_factor=(period > 20) ? period : 20;

  if (max_period < period) 
    max_period=32*period;

  // Beginning adaptive burn-in 
  while (period <= max_period)
    {

      // simulate approximately period/p_saved draws 
      Simulate(false,string(),period,reinitialize_factor,false,false);

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
	  //new_scale=(1.0 + 4.0*(previous_ratio - target_scale)/(1.0 - target_scale))*high_scale;
	  new_scale = (previous_ratio > upper_bound) ? 5.0*high_scale : (log_mid/log(previous_ratio))*high_scale;
	}
      else if (high_jump_ratio < 0.0)
	{
	  best_scale = scale; 
	  //new_scale=(0.2 + 0.8*previous_ratio/target_scale)*low_scale;
	  new_scale = (previous_ratio < lower_bound) ? 0.2*low_scale : (log_mid/log(previous_ratio))*low_scale;
	}
      else 
	{
	  best_scale = new_scale = ((target_scale-low_jump_ratio)*low_scale + (high_jump_ratio-target_scale)*high_scale)/(high_jump_ratio-low_jump_ratio);
	  period*=2; 
	  low_jump_ratio = -1.0; 
	  high_jump_ratio = -1.0; 
	}
			
      // reset scale
      scale = new_scale; 
    }	  

  // set scale
  scale=best_scale;
  
  return scale;
}

/*
   Set reinitialize_factor to -1.0 to never reinitialize.
*/
void CEquiEnergyModel::Simulate(bool if_storage, const string &sample_file_name, int number_to_save, int reinitialize_factor, bool start_from_current_sample, bool verbose)
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
  int verbose_factor = number_to_save/20;
  if (verbose_factor < 500) verbose_factor=500;
  CSampleIDWeight x_new;
  for (int verbose_count=verbose_factor, count=reinitialize_factor; nSaved < number_to_save; )
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
	  count=reinitialize_factor;
	  nRestarts++;
	}

      // write previous draw
      if (dw_uniform_rnd() <= parameter->p_save)
	{
	  nSaved++;
	  count--;
	  verbose_count--;
	  if (if_storage)
	    SaveSampleToStorage(current_sample); 
	  if (if_write_file)
	    write(output_file, &current_sample);     

	  // make equi-energy jump
	  if (dw_uniform_rnd() < parameter->pee)  
	    {
	      nEEProposed++;
	      if (storage->DrawSample(energy_level+1, current_ring, x_new))
		{
		  if (storage->BinIndex(energy_level+1,-current_sample.weight) != current_ring)
		    {
		      cerr << "DrawSample() did not draw from desired ring (ring=" << current_ring << ")" << endl;
		      abort();
		    }

		  if (log(dw_uniform_rnd()) <= (x_new.weight - current_sample.weight)*(1.0/K(energy_level) - 1.0/K(energy_level+1)))
		    {
		      nEEJumps++;
		      current_sample = x_new; 
		      current_sample.id = (int)(time(NULL)-timer_when_started);
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
	  current_sample.id = (int)(time(NULL)-timer_when_started);
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

double CEquiEnergyModel::GetMetropolisProposal(CSampleIDWeight &new_x)
{
  // int n=IndependentDirections[energy_level].cols;
  // TDenseVector v(n);
  // double density=density_constant; 
  // if (parameter->p_select >= 1.0)
  //   v=scale(energy_level)*RandomNormalVector(n);
  // else
  //   {
  //     for (int i=n-1; i >= 0; i--)
  // 	v.vector[i]=(dw_uniform_rnd() <= parameter->p_select) ? scale(energy_level) * dw_gaussian_rnd() : parameter->tiny * scale(energy_level) * dw_gaussian_rnd();     
  //   }
  // new_x.data=current_sample.data + IndependentDirections[energy_level]*v;
  // new_x.weight=target->LogPosterior(new_x.data);
  // new_x.calculated=true;

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
  stringstream convert;
  convert << parameter->run_id << "/" << parameter->run_id << ".Diagnostic." << energy_level << "." << node;
  string filename = parameter->storage_dir + convert.str();
  ofstream output_file;
  output_file.open(filename.c_str(), ios::out);
  if (!output_file)
    {
      cerr << "Error opening " << filename << endl; 
      abort(); 	
    }

  output_file << "Scale: " << scale << endl;
  output_file << "Temperature: " << K(energy_level) << endl;
  output_file << "Effective sample size: " << ess(energy_level) << endl;
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

  output_file << setprecision(4) << "Equi-energy: " << nEEJumps << "/" << nEEProposed << " = " << (nEEProposed ? 100*(double)nEEJumps/(double)nEEProposed : 0.0) << endl;
  output_file << "Metropolis-Hastings: " << nMHJumps << "/" << nMHProposed << " = " << (nMHProposed ? 100*(double)nMHJumps/(double)nMHProposed : 0.0) << endl;
  output_file << "Importance draw starts:" << nRestarts << endl;
  output_file << "Number saved: " << nSaved << endl;
  
  output_file << "Independent directions:" << endl;
  //output_file << IndependentDirections[energy_level] << endl;
    
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
  out_file << Or;
  out_file << Di;
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

  in >> K.vector[level];
  in >> scale;
  OrthonormalDirections.Resize(parameter->nParameters,parameter->nParameters);
  in >> OrthonormalDirections;
  SqrtDiagonal.Resize(parameter->nParameters);
  in >> SqrtDiagonal;
  in.close();
  iImportanceSamples=0;

  // // this provides compatibility with old code
  // IndependentDirections[level]=OrthonormalDirections;
  // for (int j=0; j < IndependentDirections[level].cols; j++)
  //   for (int i=0; i < IndependentDirections[level].cols; i++)
  //     IndependentDirections[level](i,j)*=SqrtDiagonal(j);
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
      x(j)=dw_tdistribution_rnd(nu);		
      log_density+=log(dw_tdistribution_pdf(x(j),nu)/SqrtDiagonal(j));		
      x(j)*=SqrtDiagonal(j);
    }
  x=InitialCenter + OrthonormalDirections*x;
  return log_density;
}

double CEquiEnergyModel::LogInitialDensity(const TDenseVector &x)
{
  double log_density=0.0;
  cout << OrthonormalDirections.rows << " " << OrthonormalDirections.cols << " " << x.dim << " " << InitialCenter.dim << endl; 
  TDenseVector y=TransposeMultiply(OrthonormalDirections,x-InitialCenter);
 cout << "here\n";
  for (int j=parameter->nParameters-1; j >= 0; j--)	
    log_density+=log(dw_tdistribution_pdf(y(j),nu)/SqrtDiagonal(j));		
  return log_density;
}


/*
   (1) Read draws from the top level
   (2) Compute InitialCenter, InitialOrthonormalDirections and InitialSqrtDiagonal
   (3) Returns effective sample size with respect to the next lower level

   Note: K(level-1) and nu, the degress of freedom for the t-distribution, must be
   set upon entry.
*/
double CEquiEnergyModel::AnalyzeInitialDraws(int level)
{
  int nParameters=parameter->nParameters, simulation_length=parameter->simulation_length;

  // read draws
  string filename;
  ifstream in;
  int node=0, count=0;
  TDenseVector x(nParameters), sum(nParameters,0.0);
  TDenseMatrix Osum(nParameters,nParameters,0.0);
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
	      sum+=x;
	      Osum+=OuterProduct(x,x);
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

  // set OrthonormalDirections, SqrtDiagonal, and InitialCenter
  InitialCenter=(1.0/(double)count)*sum;
  Eig(SqrtDiagonal,OrthonormalDirections,(1.0/(double)count)*Osum);
  for (int i=nParameters-1; i >= 0; i--)
    SqrtDiagonal(i)=sqrt(SqrtDiagonal(i));

  // compute ESS
  TDenseVector log_posterior(simulation_length), log_initial_density(simulation_length);
  for (int i=simulation_length-1; i >= 0; i--)
    {	
      log_initial_density(i)=InitialDraw(x);
      log_posterior(i)=target->LogPosterior(x);
    }

  return EffectiveSampleSize(K(level-1),log_posterior,log_initial_density);
}
