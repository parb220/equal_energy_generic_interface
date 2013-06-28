#include <cmath>
#include <algorithm>
#include <vector> 
#include <functional>
#include <ctime>
#include "CEqualEnergyModel.h"
#include "CIndexBlockSizeScale.h"

extern "C" {
	#include "dw_metropolis_theta.h"
	#include "dw_math.h"
	#include "dw_switchio.h"
	#include "dw_matrix.h"
	#include "dw_rand.h"
}

using namespace std;

double CEquiEnergyModel::LogPosterior(const double *x, int nx)
{
	double *old_x = new double[nx]; 
	// Save the old parameters stored in target_model to old_x
	ConvertThetaToFreeParameters(target_model, old_x); 
	ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) ); 
	
	// Post x to target_model
	ConvertFreeParametersToTheta(target_model, x); 
	ConvertFreeParametersToQ(target_model, x+NumberFreeParametersTheta(target_model) ); 
	double log_posterior = LogPosterior_StatesIntegratedOut(target_model);

	// Post old_x back to target_model 
	ConvertFreeParametersToTheta(target_model, old_x); 
	ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) ); 
	delete [] old_x; 

	return log_posterior;  
}

double CEquiEnergyModel::LogLikelihood(const double *x, int nx)
{
	double *old_x = new double(nx); 
	// Save the old parameters stored in target_model to old_x
	ConvertThetaToFreeParameters(target_model, old_x); 
	ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) ); 

	// post x to target_model
	ConvertFreeParametersToTheta(target_model, x); 
	ConvertFreeParametersToQ(target_model, x+NumberFreeParametersTheta(target_model) ); 
	double log_likelihood = LogLikelihood_StatesIntegratedOut(target_model); 

	// post old_x back to target_model 
	ConvertFreeParametersToTheta(target_model, old_x); 
	ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) ); 
	delete [] old_x; 

	return log_likelihood; 
}

int CEquiEnergyModel::EE_Draw(const CEESParameter &parameter, CStorageHead &storage)
{
	CSampleIDWeight x_new; 
	int new_sample_code = NO_JUMP; 

	if (energy_level == parameter.number_energy_level-1 || dw_uniform_rnd() > parameter.pee)
		if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data))
		{
			current_sample = x_new; 
			current_sample.id = (int)(time(NULL)-timer_when_started);
			new_sample_code = METROPOLIS_JUMP; 
		}
	else 
	{
		// draw x_new from bin of the higher level of the same energy; 
		int bin_id = parameter.BinIndex(-current_sample.weight, energy_level+1);
		if (storage.DrawSample(bin_id, x_new)) // if a sample is successfully draw from bin
		{
			// calculate log_ratio in the current and the higher levels
			double log_ratio = parameter.LogRatio_Level(-x_new.weight, -current_sample.weight, energy_level); 
			log_ratio += parameter.LogRatio_Level(-current_sample.weight(), -x_new.weight(), energy_level+1); 
			if (log(dw_uniform_rnd()) <= log_ratio)
			{
				// accept the new sample
				current_sample = x_new; 
				current_sample.id = (int)(time(NULL)-timer_when_started);
				new_sample_code = EQUI_ENERGY_JUMP; 
			}	
			else 
				if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data))
				{
					current_sample = x_new; 
					current_sample.id = (int)(time(NULL)-timer_when_started); 
					new_sample_code = METROPOLIS_JUMP; 
				}
		}
		else 
			if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data) )
			{
				current_sample = x_new; 
				current_sample.id = (int)(time(NULL)-timer_when_started); 
				new_sample_code = METROPOLIS_JUMP; 
			}
	}
	
	return new_sample_code; 
}

// A sample is randomly taken from a pool (with size desired_pool_size) of samples where the pool is formed by samples with higher log-posteriors. Note that samples with higher log-posterior values are stored in smaller-indexed bins. So, if the desired pool size is 10 while the size of the first bin is 100, then only the first bin will be used. In contrast, if the desired pool size is 100 while the total number of samples in the first 3 bins is barely greater than 100, then the first 3 bins will be used. 
bool CEquiEnergyModel::Initialize(CStorageHead &storage, unsigned int start_bin, unsigned int end_bin, size_t desiredPoolSize)
{
	size_t N=end_bin-start_bin+1; 
	unsigned int nSample_Total=0;
	vector<unsigned int>nSample_Bin(N,0);
	for (unsigned int bin=start_bin; bin<=end_bin; bin++)
	{
		nSample_Bin[bin-start_bin] = storage.GetNumberRecrod(bin); 
		nSample_Total += nSample_Bin[bin-start_bin]; 
	}

	unsigned int nLumSum = 0, bin= start_bin; 
	unsigned int random_index = dw_uniform_int(desiredPoolSize < nSample_Total ? desiredPoolSize : nSample_Total); 
	while (bin <= end_bin && !(random_index >= nLumSum && random_index < nLumSum + nSample_Bin[bin-start_bin]) )
	{	
		nLumSum += nSample_Bin[bin-start_bin]; 
		bin ++;
	}
	if (random_index >= nLumSum && random_index < nLumSum + nSample_Bin[bin-start_bin])
	{
		if(storage.DrawSample(bin, current_sample))
		{
			current_sample.id = (int)(time(NULL)-timer_when_started); 
			// Because all samples stored in storage have had their log-posterior calculated and stored 
			// together with the sample values, there is no need to recalculate log-posterior at this moment
			// current_sample.weight = LogPosterior(current_sample.data); 
                	return true; 
		}	
	}	
	return false; 	
}

bool CEquiEnergyModel::InitializeWithBestSample(CStorageHead &storage, unsigned int start_bin, unsigned int end_bin)
{
        unsigned int bin=start_bin, number_sample;
        while (bin <= end_bin && (number_sample = storage.GetNumberRecrod(bin)) <= 0)
                bin ++;
        if (bin > end_bin)
                return false;
        if (storage.DrawLeastWeightSample(bin, current_sample))
	{
		current_sample.id = (int)(time(NULL)-timer_when_started); 
                return true;
	}
	else
        	return false;
}

double CEquiEnergyModel::BurnIn(size_t burn_in_length)
{
	CSampleIDWeight x_new; 
	unsigned int nJump =0; 
	double max_posterior = current_sample.weight; 
	for (unsigned int i=0; i<burn_in_length; i++)
	{
		if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data) )
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

double CEquiEnergyModel::Simulation_Within(const CEESParameter &parameter, CStorageHead &storage, bool if_storage, FILE *dw_file_out, bool if_write_file)
{
	CSampleIDWeight x_new; 
	unsigned int nJump =0; 
	double max_posterior = current_sample.weight; 
	for (unsigned int i=0; i<parameter.simulation_length; i++)
	{
		if (i%parameter.deposit_frequency == 0)
		{
			if (if_storage)
				SaveSampleToStorage(current_sample, parameter, storage); 
			if (if_write_file)
				SaveSampleToFile(current_sample, dw_file_out); 
		
		}
		if (metropolis->BlockRandomWalkMetropolis(x_new.weight, x_new.data, current_sample.data) )
                {
                        current_sample = x_new;
                        current_sample.id = (int)(time(NULL)-timer_when_started);
                        max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior;
                        nJump ++;
                }
	}
	cout << "MH Jump " << nJump << " out of " << parameter.simulation_length << " in simulation.\n"; 
	return max_posterior;  
}

double CEquiEnergyModel::Simulation_Cross(const CEESParameter &parameter, CStorageHead &storage, bool if_storage, FILE *dw_file_out, bool if_write_file)
{
	CSampleIDWeight x_new;

	unsigned int nEEJump=0, nMHJump=0; 

	double max_posterior = current_sample.weight; 
	for (int i=0; i<parameter.simulation_length; i++)
	{
		if (i%parameter.deposit_frequency ==0)
		{
			if (if_storage)
				SaveSampleToStorage(current_sample, parameter, storage); 
			if (if_write_file)
				SaveSampleToFile(current_sample, dw_file_out); 
		}
		int jump_code = EE_Draw(parameter, storage); 
		if (jump_code == EQUI_ENERGY_JUMP)
			nEEJump++; 
		else if (jump_code == METROPOLIS_JUMP)
			nMHJump++; 
		max_posterior = max_posterior > current_sample.weight ? max_posterior : current_sample.weight; 
	}
	
	cout << "EE Jump " << nEEJump << " out of " << parameter.simulation_length << " in simulation.\n"; 
	cout << "MH Jump " << nMHJump << " out of " << parameter.simulation_length << " in simulation.\n"; 
	return max_posterior; 
}
///////////////////////////////////////////////////////////

void DetermineBlockSize(int, int *, int, int); 
void DetermineScale_ReshuffleIndex_IfNeeded(double *, int *, int , const int *, int, const vector<vector<CIndexBlockSizeScale> > &);
// random generator function:
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

void CStateModel::SaveSampleToStorage(const CSampleIDWeight &sample, const CEESParameter &parameter, CStorageHead &storage)
{
	int binIndex = parameter.BinIndex(sample.GetWeight(), energy_level);
        storage.DepositSample(binIndex, sample);
}

void CStateModel::SaveSampleToFile(const CSampleIDWeight &sample, FILE *file)
{
	if (file == NULL)
                return;
        fprintf(file, "%le %le", sample.log_prob, sample.GetWeight());  
        for (int j=0; j<thetaDim+qDim; j++)
                fprintf(file, " %le", *(sample.GetData()+j));
        fprintf(file, "\n");
}

void CStateModel::GetCurrentSample(CSampleIDWeight &sample, bool if_calculate_energy) const
{
	sample.SetDataDimension(thetaDim+qDim); 

	if (if_calculate_energy)
		LogPosterior_StatesIntegratedOut(model);

	ConvertThetaToFreeParameters(model, sample.Data());
        ConvertQToFreeParameters(model, sample.Data()+thetaDim);
	sample.SetID((int)time(NULL)); 

	if (model->if_bounded > 0)
		sample.log_prob = -(model->energy > model->h ? model->energy : model->h)/model->t;
	else
		sample.log_prob = -model->energy;

	sample.SetWeight(model->energy);
}

/*int CStateModel::Calibrate_Diag(vector<CIndexBlockSizeScale> &scale_1, vector<CIndexBlockSizeScale> &scale_n, double center_s, double center_a, int period_s, int period_a, int max_period_s, int max_period_a)
{
	TMetropolis_theta *metropolis=model->metropolis_theta;
	Setup_metropolis_theta_diagonal(metropolis);
  	Setup_metropolis_theta_blocks_single(metropolis);

	double *scale_array_1 = new double[thetaDim]; 	
	double *scale_array_n = new double[thetaDim]; 

	int ndraws=Calibrate_metropolis_theta_two_pass_scale_returned(scale_array_1, scale_array_n, thetaDim, model, center_s, period_s, max_period_s, center_a, period_a, max_period_a, 0);
	
	scale_1.resize(thetaDim); 	
	scale_n.resize(thetaDim); 
	for (int i=0; i<thetaDim; i++)
	{
		scale_1[i] = CIndexBlockSizeScale(i, 1, scale_array_1[i]); 
		scale_n[i] = CIndexBlockSizeScale(i, thetaDim, scale_array_n[i]); 
	}
	delete [] scale_array_1; 
	delete [] scale_array_n;  
	return ndraws; 
}

int CStateModel::Calibrate_Diag_Random_Block(vector<CIndexBlockSizeScale> &scale_random_block, const vector<CIndexBlockSizeScale> &scale_1, const vector<CIndexBlockSizeScale > &scale_n, double center, int period, int max_period, const CEESParameter &parameter)
// scale_random_block: scale of each dimension estimated during random-block tuning
{
	TMetropolis_theta *metropolis=model->metropolis_theta;
        Setup_metropolis_theta_diagonal(metropolis);
	
	int *index_shuffled = new int [thetaDim];
        for (int i=0; i<thetaDim; i++)
                index_shuffled [i] = i;
	random_shuffle(index_shuffled, index_shuffled + thetaDim, myrandom);

        // n_block: number of blocks
      	int n_block = ceil((double)thetaDim/(double)parameter.size_per_block);
        // block_size: number of indices in each block
        int *block_size = new int [n_block];
        DetermineBlockSize(thetaDim, block_size, n_block, parameter.size_per_block);
		
	double *scale_array_random_block = new double[thetaDim]; 
	double *scale_array_1 = new double[thetaDim]; 
	double *scale_array_n = new double[thetaDim]; 
	for (int k=0; k<thetaDim; k++)
	{
		scale_array_1[k] = scale_1[k].GetScale(); 
		scale_array_n[k] = scale_n[k].GetScale(); 
	} 
	int ndraws = Calibrate_metropolis_theta_random_block(scale_array_random_block, scale_array_1, scale_array_n, thetaDim, index_shuffled, thetaDim, block_size, n_block, model, center, period, max_period, 1); 
	
	scale_random_block.resize(thetaDim); 
	int block_offset = 0; 
	for (int i_block=0; i_block<n_block; i_block++)
	{
		for (int k=block_offset; k<block_offset+block_size[i_block]; k++)
                {
                        if (block_size[i_block] == 1)
                                scale_random_block[index_shuffled[k]] = scale_1[index_shuffled[k]];
                        else if (block_size[i_block] == thetaDim)
                                scale_random_block[index_shuffled[k]] = scale_n[index_shuffled[k]];
                        else
                                scale_random_block[index_shuffled[k]] = CIndexBlockSizeScale(index_shuffled[k], block_size[i_block], scale_array_random_block[index_shuffled[k]]);
                }
		block_offset += block_size[i_block]; 
	}

	delete []scale_array_random_block;
	delete []scale_array_1;
	delete []scale_array_n;   
	delete []index_shuffled; 
	delete []block_size; 
	return ndraws; 
}

int CStateModel::Calibrate_Variance_Initialize_Sample_Random_Block(vector<CIndexBlockSizeScale> &scale_random_block, const vector <CIndexBlockSizeScale> &scale_1, const vector<CIndexBlockSizeScale> &scale_n, string sample_file, double center, int period, int max_period, const CEESParameter &parameter)
// scale_random_block: scale of each dimension estimated during random-block tuning
{
	TMatrix posterior_draws = (TMatrix)NULL;
	posterior_draws = dw_ReadPosteriorDraws((FILE *)NULL, (char*)sample_file.c_str(), (char *)NULL, thetaDim+qDim); 
	if (!posterior_draws)
	{
		scale_random_block.clear(); 
		return -1; // fail
	}
	
	TMatrix variance = (TMatrix)NULL, X=(TMatrix)NULL;  
	TVector mean = (TVector)NULL, x=(TVector)NULL, theta=(TVector)NULL; 

	// Read out previous draws and calculate variance
	InitializeMatrix(variance = CreateMatrix(thetaDim, thetaDim), 0.0); 
	X=CreateMatrix(thetaDim, thetaDim); 
	InitializeVector(mean = CreateVector(thetaDim), 0.0); 
	x=CreateVector(thetaDim+qDim+2); 
	theta=CreateVector(thetaDim); 
	for (int i=RowM(posterior_draws)-1; i>=0; i--)
	{
		RowVector(x, posterior_draws, i); 
		SubVector(theta, x, 2, thetaDim); 
		AddVV(mean, mean, theta); 
		OuterProduct(X, theta, theta); 
		AddMM(variance, variance, X); 
	}
	ProductVS(mean, mean, 1.0/(double)RowM(posterior_draws)); 
	ProductMS(variance, variance, 1.0/(double)RowM(posterior_draws)); 
	SubtractMM(variance, variance, OuterProduct(X, mean, mean)); 
	
	// Metropolis calibration
	TMetropolis_theta *metropolis=model->metropolis_theta;
        Setup_metropolis_theta_variance(metropolis, variance);

	int *index_shuffled = new int [thetaDim];
        for (int i=0; i<thetaDim; i++)
                index_shuffled [i] = i;
	random_shuffle(index_shuffled, index_shuffled + thetaDim, myrandom);
	// n_block: number of blocks
	int n_block = ceil((double)thetaDim/(double)parameter.size_per_block);
	// block_size: number of indices in each block
	int *block_size = new int [n_block];
	DetermineBlockSize(thetaDim, block_size, n_block, parameter.size_per_block);
	
	double *scale_array_random_block = new double [thetaDim];
	double *scale_array_1 = new double[thetaDim]; 
	double *scale_array_n = new double[thetaDim];  
        for (int k=0; k<thetaDim; k++)
	{
		scale_array_1[k] = scale_1[k].GetScale(); 
		scale_array_n[k] = scale_n[k].GetScale(); 
	}
	int ndraws=Calibrate_metropolis_theta_random_block(scale_array_random_block, scale_array_1, scale_array_n, thetaDim, index_shuffled, thetaDim, block_size, n_block, model, center, period, max_period, 1);
	
	scale_random_block.resize(thetaDim); 
	int block_offset = 0;
	for (int i_block=0; i_block<n_block; i_block++)
	{
		for (int k=block_offset; k<block_offset+block_size[i_block]; k++)
		{
			if (block_size[i_block] == 1)
				scale_random_block[index_shuffled[k]] = scale_1[index_shuffled[k]]; 
			else if (block_size[i_block] == thetaDim)
				scale_random_block[index_shuffled[k]] = scale_n[index_shuffled[k]]; 
			else 
				scale_random_block[index_shuffled[k]] = CIndexBlockSizeScale(index_shuffled[k], block_size[i_block], scale_array_random_block[index_shuffled[k]]);
		}
		block_offset += block_size[i_block];
	}

	delete []scale_array_random_block; 
	delete []scale_array_1; 
	delete []scale_array_n;
	delete []index_shuffled; 
	delete []block_size; 
	// free space
	FreeMatrix(X); 
	FreeMatrix(variance); 
	FreeVector(x); 
	FreeVector(mean); 
	FreeVector(theta); 
        FreeMatrix(posterior_draws); 
	return ndraws;
}

int CStateModel::Calibrate_Variance_Initialize_Sample(vector<CIndexBlockSizeScale> &scale_1, vector<CIndexBlockSizeScale> &scale_n, string sample_file, double center_s, double center_a, int period_s, int period_a, int max_period_s, int max_period_a)
{
	TMatrix posterior_draws = (TMatrix)NULL;
	posterior_draws = dw_ReadPosteriorDraws((FILE *)NULL, (char*)sample_file.c_str(), (char *)NULL, thetaDim+qDim); 
	if (!posterior_draws)
	{
		scale_1.clear(); 
		scale_n.clear(); 
		return -1; // fail
	}
	
	TMatrix variance = (TMatrix)NULL, X=(TMatrix)NULL;  
	TVector mean = (TVector)NULL, x=(TVector)NULL, theta=(TVector)NULL; 

	// Read out previous draws and calculate variance
	InitializeMatrix(variance = CreateMatrix(thetaDim, thetaDim), 0.0); 
	X=CreateMatrix(thetaDim, thetaDim); 
	InitializeVector(mean = CreateVector(thetaDim), 0.0); 
	x=CreateVector(thetaDim+qDim+2); 
	theta=CreateVector(thetaDim); 
	for (int i=RowM(posterior_draws)-1; i>=0; i--)
	{
		RowVector(x, posterior_draws, i); 
		SubVector(theta, x, 2, thetaDim); 
		AddVV(mean, mean, theta); 
		OuterProduct(X, theta, theta); 
		AddMM(variance, variance, X); 
	}
	ProductVS(mean, mean, 1.0/(double)RowM(posterior_draws)); 
	ProductMS(variance, variance, 1.0/(double)RowM(posterior_draws)); 
	SubtractMM(variance, variance, OuterProduct(X, mean, mean)); 
	
	// Metropolis calibration
	TMetropolis_theta *metropolis=model->metropolis_theta;
        Setup_metropolis_theta_variance(metropolis, variance);
        Setup_metropolis_theta_blocks_single(metropolis);

	double *scale_array_1 = new double [thetaDim]; 
	double *scale_array_n = new double [thetaDim]; 
        int ndraws=Calibrate_metropolis_theta_two_pass_scale_returned(scale_array_1, scale_array_n, thetaDim, model, center_s, period_s, max_period_s, center_a, period_a, max_period_a, 0);
	
	scale_1.resize(thetaDim); 
	scale_n.resize(thetaDim); 
	for (int i=0; i<thetaDim; i++)
	{
		scale_1[i] = CIndexBlockSizeScale(i, 1, scale_array_1[i]); 
		scale_n[i] = CIndexBlockSizeScale(i, thetaDim, scale_array_n[i]); 
	}
	
	delete []scale_array_1;
	delete []scale_array_n;  
	// free space
	FreeMatrix(X); 
	FreeMatrix(variance); 
	FreeVector(x); 
	FreeVector(mean); 
	FreeVector(theta); 
        FreeMatrix(posterior_draws); 
	return ndraws;
}*/

bool CStateModel::LoadMetropolisParameter(string theta_file)
{
	TMetropolis_theta *metropolis=model->metropolis_theta;
	if (Read_metropolis_theta((FILE*)NULL, (char*)theta_file.c_str(), (char*)NULL, metropolis) )
		return true; 
	else 
		return false; 
}


void DetermineScale_ReshuffleIndex_IfNeeded(double *scale_array, int *index_shuffled, int n_array_index, const int *block_size, int n_block, const vector<vector<CIndexBlockSizeScale> > &scale)
{
	bool continue_reshuffle = true, continue_check; 
	int block_offset; 
	while(continue_reshuffle)
	{
		random_shuffle(index_shuffled, index_shuffled+n_array_index, myrandom);
		continue_reshuffle = false;
		continue_check = true;  
		// check if scale has been specified for a particular dimension and block_size
		// if not, then needs to reshuffle again
		block_offset = 0; 
		int i_block = 0;  
		while (continue_check && i_block<n_block)
		{
			int i_dim = block_offset; 
			while (continue_check && i_dim<block_offset+block_size[i_block])
			{
				if( scale[index_shuffled[i_dim]][block_size[i_block]-1].GetScale()<0)
				{
					continue_check = false; 
					continue_reshuffle=true; 
				}
				else 
				{
					if (block_size[i_block] == 1 || block_size[i_block] == n_array_index) 
						scale_array[index_shuffled[i_dim]] = scale[index_shuffled[i_dim]][block_size[i_block]-1].GetScale(); 
					else 
						scale_array[index_shuffled[i_dim]] = scale[index_shuffled[i_dim]][block_size[i_block]-1].GetScale(); // > scale[index_shuffled[i_dim]][n_array_index-1].GetScale() ? scale[index_shuffled[i_dim]][block_size[i_block]-1].GetScale() *0.1 : scale[index_shuffled[i_dim]][n_array_index-1].GetScale(); 
				}
				i_dim++; 
			}
			block_offset += block_size[i_block];
			i_block++; 
		}
	}
}

double CStateModel::BurnIn_MultipleBlock(const vector<vector<CIndexBlockSizeScale> > &scale, const CEESParameter &parameter, int burn_in_length)
{
	// index_shuffled: original 0 ... (thetaDim-1), then randomly shuffeled once a while 
	int *index_shuffled = new int [thetaDim];
        for (int i=0; i<thetaDim; i++)
                index_shuffled [i] = i;

	// n_block: number of blocks
	int n_block = ceil((double)thetaDim/(double)parameter.size_per_block);

	// block_size: number of indices in each block
	int *block_size = new int [n_block];
	DetermineBlockSize(thetaDim, block_size, n_block, parameter.size_per_block); 

	// scales to be used
	double *scale_array = new double[thetaDim]; 
	DetermineScale_ReshuffleIndex_IfNeeded(scale_array, index_shuffled, thetaDim, block_size, n_block, scale); 

	int nJump = 0; 
	nJump += DrawAllMultipleBlock(scale_array, index_shuffled, block_size, n_block, model);  
	double min_energy = model->energy; 
	for (int i=0; i<burn_in_length; i++)
	{
		if (i%(parameter.shuffle_frequency * parameter.deposit_frequency) == 0)
			// random_shuffle(index_shuffled, index_shuffled + thetaDim, myrandom); 
			DetermineScale_ReshuffleIndex_IfNeeded(scale_array, index_shuffled, thetaDim, block_size, n_block, scale); 
		

		nJump+= DrawAllMultipleBlock(scale_array, index_shuffled, block_size, n_block, model); 
		min_energy = min_energy <= model->energy ? min_energy : model->energy; 
	}	
	cout << "MH Jump " << nJump << " out of " << burn_in_length << " in burning-in.\n"; 
	delete [] scale_array;
	delete [] index_shuffled; 
	delete [] block_size; 
	return min_energy; 
}


void *CStateModel::objfun(int *mode, int *n, double *x, double *f, double *g, int *state)
{
	//* paste x to model
	ConvertFreeParametersToTheta(model, x); 
	ConvertFreeParametersToQ(model, x+thetaDim); 
	LogPosterior_StatesIntegratedOut(model); 
	// Since npsol is to minimize the objective function, 
	// we here use energy (instead of log_prob) as the 
	// objective function
	if (model->if_bounded > 0)
		*f = (model->energy > model->h ? model->energy : model->h)/model->t; 
	else 
		*f = model->energy; //*/
}

void *CStateModel::confun(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate)
{ 
	*mode = -1;
}

double CStateModel::objfun_Csminwel(double *x, int n, double **args, int *dims)
{
	//*
	ConvertFreeParametersToTheta(model, x);
	ConvertFreeParametersToQ(model, x+thetaDim);
	LogPosterior_StatesIntegratedOut(model);
	
	/*double value; 
	if (model->if_bounded > 0)
		value = (model->energy > model->h ? model->energy : model->h)/model->t;	
	else 
		value = model->energy; 
	return value; 
	*/
	return model->energy; 
	
	/*double value=0.0; 
	for (int i=0; i<thetaDim+qDim; i++)
		value += (x[i]-i)*(x[i]-i); 
	return value; */
}

void DetermineBlockSize(int n, int *block_size, int n_block, int size_per_block)
{
	if ( n%size_per_block)
	{
                block_size[0] = n%size_per_block;
                for (int k=1; k<n_block; k++)
                        block_size[k] = size_per_block;
        }
	else
        {
                for (int k=0; k<n_block; k++)
                        block_size[k] = size_per_block;
        }
}

double CStateModel::Simulation_Within_MultipleBlock(const vector<vector<CIndexBlockSizeScale> >&scale, const CEESParameter &parameter, CStorageHead &storage, bool if_storage, FILE *dw_file_out, bool if_write_file, const gsl_rng *r)
{
	// index_shuffled: original 0 ... (thetaDim-1), then randomly shuffeled once a while 
	int *index_shuffled = new int [thetaDim];
        for (int i=0; i<thetaDim; i++)
                index_shuffled [i] = i;

	// n_block: number of blocks
	int n_block = ceil((double)thetaDim/(double)parameter.size_per_block);

	// block_size: number of indices in each block
	int *block_size = new int [n_block];
	DetermineBlockSize(thetaDim, block_size, n_block, parameter.size_per_block); 
	
	// scale: scale for each dimension and a particular block size
	double *scale_array = new double[thetaDim]; 
	DetermineScale_ReshuffleIndex_IfNeeded(scale_array, index_shuffled, thetaDim, block_size, n_block, scale); 
	
	// Simulation started
	CSampleIDWeight sample; 
	int nJump = 0; 

        nJump += DrawAllMultipleBlock(scale_array, index_shuffled, block_size, n_block, model);
	double min_energy = model->energy;
	for (int i=0; i<parameter.simulation_length; i++)
        {
                if (i%parameter.deposit_frequency == 0)
                {
                        GetCurrentSample(sample, false); // no energy calculation
			if (if_storage)
                       		SaveSampleToStorage(sample, parameter, storage);
                        if (if_write_file)
				SaveSampleToFile(sample, dw_file_out);
                }
                if (i%(parameter.shuffle_frequency * parameter.deposit_frequency) == 0)
                        // random_shuffle(index_shuffled, index_shuffled+thetaDim, myrandom);
                	DetermineScale_ReshuffleIndex_IfNeeded(scale_array, index_shuffled, thetaDim, block_size, n_block, scale);


		// DetermineBlockSize(thetaDim, block_size, n_block, parameter.size_per_block); 
		// It is very surprising that block_size is changed during the call of DrawAll
		// MultipleBlock (DrawThetaMultipleBlock, Draw_theta_multiple_block). So we have
		// to recall DetermineBlockSize to have block_size reset to the desired values
                nJump += DrawAllMultipleBlock(scale_array, index_shuffled, block_size, n_block, model);
                min_energy = min_energy <= model->energy ? min_energy : model->energy;
        }

	cout << "MH Jump " << nJump << " out of " << parameter.simulation_length << " in simulation.\n"; 	

	delete [] scale_array;
	delete [] index_shuffled; 
	delete [] block_size; 
	return min_energy; 
}

int CStateModel::EE_Draw_MultipleBlock(const double *scale_array, int *shuffled_index, int *block_size, int n_block, const CEESParameter &parameter, CStorageHead &storage, const gsl_rng *r)
{
	int jump_code = -1; 

	double uniform_draw = gsl_rng_uniform(r); 
	CSampleIDWeight x_current, x_new;  

	if (energy_level == parameter.number_energy_level-1 || uniform_draw > parameter.pee)
		jump_code = DrawAllMultipleBlock(scale_array, shuffled_index, block_size, n_block, model); 
	else 
	{
		// The current state is either from DrawAll (or DrawAllMultipleBlock) or through equi-energy sampling. 
		// If it is the former case, then energy is already calculated inside the call. If it is the latter case, 
		// then the current state is the same as a previously drawn sample, whose energy is also calculated already.
		// Therefore, there is no need to calculate energy here.
		GetCurrentSample(x_current, false); 
		
		// draw x_new from bin of the higher temp level of the same energy
		int bin_id = parameter.BinIndex(x_current.GetWeight(), energy_level+1); 
		if (storage.DrawSample(bin_id, r, x_new))
		{
			// log-ratio in the current and the higher levels.
			double log_ratio = parameter.LogRatio_Level(x_new.GetWeight(), x_current.GetWeight(), energy_level); 
			log_ratio += parameter.LogRatio_Level(x_current.GetWeight(), x_new.GetWeight(), energy_level+1); 
			uniform_draw = gsl_rng_uniform(r); 
			if (log(uniform_draw) <= log_ratio)
			{
				// accept x_new as x_current, however the log_prob needs to be recalculated based on current level's threshold
				x_current = x_new; 
				x_current.log_prob = -(x_current.GetWeight() > parameter.h[energy_level] ? x_current.GetWeight() : parameter.h[energy_level])/parameter.t[energy_level]; 
				SetAsCurrentSample(x_current, false); 
				jump_code = -1; 
			}
			else 
				jump_code = DrawAllMultipleBlock(scale_array, shuffled_index, block_size, n_block, model); 
		}
		else 
			jump_code = DrawAllMultipleBlock(scale_array, shuffled_index, block_size, n_block, model); 
	}
	return jump_code;  
}

double CStateModel::Simulation_Cross_MultipleBlock(const vector<vector<CIndexBlockSizeScale> > &scale, const CEESParameter &parameter, CStorageHead &storage, bool if_storage, FILE *dw_file_out, bool if_write_file, const gsl_rng *r)
{
	// scale_1: scales when each dimension corresponds to a block
	// scale_n: scales when all dimensions correspond to a block

	// index_shuffled: contains randomly shuffled dimension indices
	int nMHJump =0, nEEJump = 0; 
	int *index_shuffled = new int [thetaDim];
        for (int i=0; i<thetaDim; i++)
                index_shuffled [i] = i;

	// n_block: number of blocks
	int n_block = ceil((double)thetaDim/(double)parameter.size_per_block); 

	// block_size: number of dimension indices in each block
	int *block_size = new int [n_block];
	DetermineBlockSize(thetaDim, block_size, n_block, parameter.size_per_block);	
	
	// scale
	double *scale_array = new double[thetaDim]; 
	DetermineScale_ReshuffleIndex_IfNeeded(scale_array, index_shuffled, thetaDim, block_size, n_block, scale);

	// Simulation started
	CSampleIDWeight sample; 
	double min_energy; 

	int jump_code = EE_Draw_MultipleBlock(scale_array, index_shuffled, block_size, n_block, parameter, storage, r);
	if (jump_code < 0)
		nEEJump ++; 
	else 
		nMHJump += jump_code; 
        min_energy = model->energy;
        for (int i=0; i<parameter.simulation_length; i++)
        {
                if (i%parameter.deposit_frequency == 0)
                {
                        GetCurrentSample(sample, false); // no energy calculation
			if (if_storage)
                        	SaveSampleToStorage(sample, parameter, storage);
			if (if_write_file)
                       		SaveSampleToFile(sample, dw_file_out);

                        if (i%(parameter.shuffle_frequency * parameter.deposit_frequency) == 0)
                                random_shuffle(index_shuffled, index_shuffled+thetaDim, myrandom);
                }

                jump_code = EE_Draw_MultipleBlock(scale_array, index_shuffled, block_size, n_block, parameter, storage, r);
		if (jump_code < 0)
			nEEJump++; 
		else 
			nMHJump += jump_code; 
                min_energy = min_energy <= model->energy ? min_energy : model->energy;

        }

	cout << "EE Jump " << nEEJump << " out of " << parameter.simulation_length << " in simulation.\n"; 
	cout << "MH Jump " << nMHJump << " out of " << parameter.simulation_length << " in simulaiton.\n"; 
	delete [] scale_array; 
	delete [] index_shuffled; 
	delete [] block_size; 
}
