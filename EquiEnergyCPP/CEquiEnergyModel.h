#ifndef _CLASS_STATE_MODEL
#define _CLASS_STATE_MODEL

#include <vector>
#include <ctime>
#include <cstdio>
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"

using namespace std;

const int METROPOLIS_JUMP = 1; 
const int EQUI_ENERGY_JUMP = -1; 
const int NO_JUMP = 0; 

class CSampleIDWeight;
class CStorageHead;
class CEESParameter;
class CMetropolis; 
class CEquiEnergyModel; 
class MinusLogPosterior_NPSOL;
class TDenseVector;
class TDenseMatrix;

// class MinusLogPosterior_CSMINWEL;

class MinusLogPosterior_NPSOL
{
public:
        static CEquiEnergyModel *model;
        static void *function(int *mode, int *n, double *x, double *f, double *g, int *nstate);
};

/*class MinusLogPosterior_CSMINWEL
{
public:
        static CEquiEnergyModel *model;
        static double function(double *x, int n, double **args, int *dims);
};*/

class CEquiEnergyModel
{
private:
	CEquiEnergyModel(const CEquiEnergyModel &);
	std::vector<TDenseVector> gmm_mean; 
	std::vector<TDenseMatrix> gmm_covariance_sqrt; 
	std::vector<double> gmm_covariance_sqrt_log_determinant; 
	std::vector<TDenseMatrix> gmm_covariance_sqrt_inverse; 

public:
///////////////////////////////////////////////////////////////
// Parameters for equi-energy sampling
	bool if_bounded; 
 	int energy_stage;		
	double t_bound; 

//////////////////////////////////////////////////////////////
// Current sample holding the following 
// 	Sample: 
// 	ID:	unique id of the sample (by default is the difference
// 		between when it is drawn and when the program is started)
// 	Weight:	energy (log posterior)
	CSampleIDWeight current_sample; 
	int timer_when_started; 
	
public:
	CMetropolis *metropolis; 	// pointer to CMetropolis
	CEESParameter *parameter;	// pointer to CEESParameter 
	CStorageHead *storage; 		// pointer to CStorageHead
	virtual double log_posterior_function(CSampleIDWeight &x)=0;
        virtual double log_likelihood_function(const CSampleIDWeight &x)=0;
	virtual double log_prior_function(const CSampleIDWeight &x) = 0; 
	virtual double log_posterior_function(const double *x, int n)=0; 
        virtual double log_likelihood_function(const double *x, int n)=0;
	virtual double log_prior_function(const double *x, int n) = 0;  

	// Draw samples
	int EE_Draw(int MH_thin); 	// equi-energy draw
public:
	double BurnIn(int burn_in_length);		// returns the maximum posteior during burn-in
	bool Initialize(int desired_pool_size, int stage);	// Initialize model (setting values for current_sample) using bins indexed from start_bin through (including) end_bin. 
	bool InitializeWithBestSample(int stage); 		// Initialize model using the best sample in the bins indexed from start_bin to end_bin
	bool InitializeWith_Kth_BestSample(int K, int stage_index);
	bool Initialize_RandomlyPickFrom_K_BestSample(int K, int stage_index); 
	bool Initialize_KMeansClustering(int K, int stage_index, vector<CSampleIDWeight> &centers) const; 
	bool Initialize_MostDistant_WithinPercentile(int K, int stage_index, vector<CSampleIDWeight > &starters, double percentile=0.50) const; 
	bool Initialize_MostDistant_WithinPercentileBand(int K, int stage_index, vector<CSampleIDWeight > &starters, double percentile=0.50) const; 
	bool Initialize_WeightedSampling(int K, int stage_index, vector<CSampleIDWeight> &starters); 
	bool InitializeFromFile(const string &file_name); 

	void Simulation_Within(bool if_storage, const string &sample_file_name=string()); 	// Simulation within the same energy stage (no jumping across stages). Returns the maximum posterior during simulation
	void Simulation_Cross(bool if_storage, const string &sample_file_name=string()); 	// Simulation across stages. Returns the maximum posterior during simulation 	
	void Simulation_Prior(bool if_storage, const string &sample_file_name); 		// Simulation from prior
	// Reweight samples
	vector<double> Reweight(const vector<CSampleIDWeight> &samples, int current_stage, int previous_stage); 

//////////////////////////////////////////////////////////////////////////////////////////
// 	make a equi-energy jump
	bool MakeEquiEnergyJump(CSampleIDWeight &y_end, const CSampleIDWeight &y_initial); 

///////////////////////////////////////////////////////////////////////////////////////////
// IO 
public:
	virtual void SaveSampleToStorage(const CSampleIDWeight &sample); 
	virtual void Take_Sample_Just_Drawn_From_Storage(const CSampleIDWeight &x); 

///////////////////////////////////////////////////////////////////////////////////////////
// Construction & Destruction functions here
public:
	CEquiEnergyModel(); 
	CEquiEnergyModel(bool _if_bounded, int eL, double _t, const CSampleIDWeight & _x=CSampleIDWeight(), time_t _time=time(NULL), CMetropolis *_metropolis =NULL, CEESParameter *_parameter=NULL, CStorageHead *_storage = NULL); 
	virtual ~CEquiEnergyModel() {} 

//////////////////////////////////////////////////////////////////////////////////////////
// HillClimb
public: 
	virtual double HillClimb_NPSOL(int nSolution, int =10, int =10, double = 1.0, double = 1000.0); 
	// virtual double HillClimb_CSMINWEL(int nSolution); 
	virtual bool DrawParametersFromPrior(double *x) const = 0; 
	double GMM_Simulation(int simulation_length);
	bool StudentT_DrawSample(CSampleIDWeight &y); 
	double StudentT_LogPDF(const CSampleIDWeight &x) const; 
	double StudentT_LogRatio(const CSampleIDWeight &x, const CSampleIDWeight &y) const; 
	bool Cauchy_DrawSample(CSampleIDWeight &y); 
	double Cauchy_LogPDF(const CSampleIDWeight &x) const; 
	double Cauchy_LogRatio(const CSampleIDWeight &x, const CSampleIDWeight &y) const; 
	bool GMM_DrawSample(CSampleIDWeight &y); 
	double GMM_LogPDF(const CSampleIDWeight &x) const; 
	double GMM_LogRatio(const CSampleIDWeight &x, const CSampleIDWeight &y) const; 
	bool WriteGaussianMixtureModelParameters(const string &file_name) const; 
	bool ReadGaussianMixtureModelParameters(const string &file_name); 
	bool AggregateGaussianMixtureModelParameters(const string &file_name);
	void ClearGaussianMixtureModelParameters(); 
	void KeepOptimalGaussianMixtureModelParameters(); 
	bool ConsolidateSampleForCovarianceEstimation(const vector<string> &filenames_merge, const string &variance_file); 
	const TDenseVector & GetGMM_Mean(int i) const { return gmm_mean[i]; }
	const TDenseMatrix & GetGMM_CovarianceSqrt(int i) const { return gmm_covariance_sqrt[i]; }
	double GetGMM_CovarianceSqrtLogDeterminant(int i) const { return gmm_covariance_sqrt_log_determinant[i]; }
	const TDenseMatrix & GetGMM_CovarianceSqrtInverse(int i) const { return gmm_covariance_sqrt_inverse[i]; }
	
friend class MinusLogPosterior_NPSOL; 
// friend class MinusLogPosterior_CSMINWEL; 
};

#endif
