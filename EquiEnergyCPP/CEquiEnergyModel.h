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


class MinusLogPosterior_NPSOL
{
public:
        static CEquiEnergyModel *model;
        static void *function(int *mode, int *n, double *x, double *f, double *g, int *nstate);
};


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
	double lambda; 

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
	int EE_Draw(); 	// equi-energy draw
	bool MakeEquiEnergyJump(CSampleIDWeight &y_end, const CSampleIDWeight &y_initial); 
	virtual void SaveSampleToStorage(const CSampleIDWeight &sample); 
	virtual void Take_New_Sample_As_Current_Sample(const CSampleIDWeight &x); 

public:
	std::vector<int> BurnIn(int burn_in_length);		// returns the maximum posteior during burn-in
	std::vector<int>Simulation_Within(bool if_storage, const string &sample_file_name=string()); 	// Simulation within the same energy stage (no jumping across stages). Returns the maximum posterior during simulation
	std::vector<int> Simulation_Cross(bool if_storage, const string &sample_file_name=string()); 	// Simulation across stages. Returns the maximum posterior during simulation 	
	std::vector<int> Simulation_Prior(bool if_storage, const string &sample_file_name); 		// Simulation from prior
	virtual bool DrawParametersFromPrior(double *x) const = 0; 
	
	std::vector<CSampleIDWeight> Initialize_MostDistant_WithinPercentile(const std::vector<CSampleIDWeight> &, int K, int stage_index, double percentile=0.50) const; 
	std::vector<CSampleIDWeight> Initialize_WeightedSampling(const std::vector<CSampleIDWeight> &, int K, int stage_index) const; 

	// Reweight samples
	vector<double> Reweight(const vector<CSampleIDWeight> &samples, int current_stage, int previous_stage) const; 


///////////////////////////////////////////////////////////////////////////////////////////
// Construction & Destruction functions here
public:
	CEquiEnergyModel(); 
	CEquiEnergyModel(bool _if_bounded, int eL, double _t, const CSampleIDWeight & _x=CSampleIDWeight(), time_t _time=time(NULL), CMetropolis *_metropolis =NULL, CEESParameter *_parameter=NULL, CStorageHead *_storage = NULL); 
	virtual ~CEquiEnergyModel() {} 

//////////////////////////////////////////////////////////////////////////////////////////
// HillClimb
public: 
	virtual double HillClimb_NPSOL(int nSolution, int =10, int =10, double = 1.0, double = 1000.0, const TDenseVector &start_point=TDenseVector(0)); 
	bool WriteGaussianMixtureModelParameters(const string &file_name) const; 
	bool WriteGaussianMixtureModelMeanAscii(const string &file_name) ; 
	bool ReadGaussianMixtureModelParameters(const string &file_name); 
	bool AggregateGaussianMixtureModelParameters(const string &file_name);
	void ClearGaussianMixtureModelParameters(); 
	void KeepOptimalGaussianMixtureModelParameters(); 
	const TDenseVector & GetGMM_Mean(int i) const { return gmm_mean[i]; }
	const TDenseMatrix & GetGMM_CovarianceSqrt(int i) const { return gmm_covariance_sqrt[i]; }
	double GetGMM_CovarianceSqrtLogDeterminant(int i) const { return gmm_covariance_sqrt_log_determinant[i]; }
	const TDenseMatrix & GetGMM_CovarianceSqrtInverse(int i) const { return gmm_covariance_sqrt_inverse[i]; }
	
friend class MinusLogPosterior_NPSOL; 
};

#endif
