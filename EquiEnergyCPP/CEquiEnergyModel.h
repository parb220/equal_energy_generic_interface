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

public:
///////////////////////////////////////////////////////////////
// Parameters for equi-energy sampling
// 	energy_level: i
// 	energy_bound: h[i]
// 	temperature_bound: t[i]
	bool if_bounded; 
 	int energy_level;		
	double t_bound; 

//////////////////////////////////////////////////////////////
// Current sample holding the following 
// 	Sample: 
// 	ID:	unique id of the sample (by default is the difference
// 		between when it is drawn and when the program is started)
// 	Weight:	energy (log posterior)
	CSampleIDWeight current_sample; 
	time_t timer_when_started; 
	
public:
	CMetropolis *metropolis; 	// pointer to CMetropolis
	CEESParameter *parameter;	// pointer to CEESParameter 
	CStorageHead *storage; 		// pointer to CStorageHead
	virtual double log_posterior_function(CSampleIDWeight &x)=0;
        virtual double log_likelihood_function(const CSampleIDWeight &x)=0;
	virtual double log_posterior_function(const double *x, int n)=0; 
        virtual double log_likelihood_function(const double *x, int n)=0; 

	// Draw samples
	int EE_Draw(int MH_thin); 	// equi-energy draw
public:
	double BurnIn(int burn_in_length);		// returns the maximum posteior during burn-in
	bool Initialize(int desired_pool_size, int level);	// Initialize model (setting values for current_sample) using bins indexed from start_bin through (including) end_bin. 
	bool InitializeWithBestSample(int level); 		// Initialize model using the best sample in the bins indexed from start_bin to end_bin
	bool InitializeWith_Kth_BestSample(int K, int level_index);
	bool Initialize_RandomlyPickFrom_K_BestSample(int K, int level_index); 
	bool Initialize_KMeansClustering(int K, int level_index, vector<CSampleIDWeight> &centers) const; 
	bool Initialize_MostDistant_WithinPercentile(int K, int level_index, vector<CSampleIDWeight > &starters, double percentile=0.50) const; 
	bool Initialize_MostDistant_WithinPercentileBand(int K, int level_index, vector<CSampleIDWeight > &starters, double percentile=0.50) const; 
	bool InitializeFromFile(const string &file_name); 

	double Simulation_Within(bool if_storage, const string &sample_file_name=string()); 	// Simulation within the same energy level (no jumping across levels). Returns the maximum posterior during simulation
	double Simulation_Cross(bool if_storage, const string &sample_file_name=string()); 	// Simulation across levels. Returns the maximum posterior during simulation 	

//////////////////////////////////////////////////////////////////////////////////////////
// 	make a equi-energy jump
	bool MakeEquiEnergyJump(CSampleIDWeight &y_end, const CSampleIDWeight &y_initial) const; 

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
	virtual double HillClimb_NPSOL(int nSolution, std::vector<TDenseVector> &, std::vector<TDenseMatrix> &, int =10, int =10, double = 1.0); 
	// virtual double HillClimb_CSMINWEL(int nSolution); 
	virtual bool DrawParametersFromPrior(double *x) const = 0; 
	double GMM_Simulation(int simulation_length, const std::vector<TDenseVector> &gm_mean, const std::vector<TDenseMatrix> &gm_covariance_sqrt);

friend class MinusLogPosterior_NPSOL; 
// friend class MinusLogPosterior_CSMINWEL; 
};

#endif
