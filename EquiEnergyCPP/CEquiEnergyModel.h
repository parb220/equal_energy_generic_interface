#ifndef _CLASS_STATE_MODEL
#define _CLASS_STATE_MODEL

#include <vector>
#include <ctime>
#include <cstdio>
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CMetropolis.h"

extern "C" {
	#include "dw_switch.h"
}

using namespace std;

const int METROPOLIS_JUMP = 1; 
const int EQUI_ENERGY_JUMP = -1; 
const int NO_JUMP = 0; 

class CMetropolis; 
class MinusLogPosterior_NPSOL; 

class CEquiEnergyModel
{
public:
///////////////////////////////////////////////////////////////
// Parameters for equi-energy sampling
// 	energy_level: i
// 	energy_bound: h[i]
// 	temperature_bound: t[i]
	bool if_bounded; 
 	unsigned int energy_level;		 
	double h_bound, t_bound;
//////////////////////////////////////////////////////////////
// Current sample holding the following 
// 	Sample: 
// 	ID:	unique id of the sample (by default is the difference
// 		between when it is drawn and when the program is started)
// 	Weight:	energy (log posterior)
	CSampleIDWeight current_sample; 

	TStateModel *target_model; 	// pointer to the TStateModel 
	CMetropolis *metropolis;	// pointer to metropolis

protected:
	time_t timer_when_started; 

///////////////////////////////////////////////////////////////
// Calculate log-posterior and log-likelihood of a given sample
// Depending on TStateModel, LogPosterior and LogLikelihood can be implemented 
// differently. 
protected:
	virtual double log_posterior_function(const double *x, size_t n); 
	virtual double log_likelihood_function(const double *x, size_t n); 
public:
	virtual double log_posterior_function(CSampleIDWeight &x); 
	// x cannot be constant because x.weight will be set as the real log_posterior calculated
	// from target_model, where the returning value is the bounded log_posterior
	// that is, return value = -(-x.weight>h_bound ? -x.weight:h_bound)/t_bound; 
	virtual double log_likelihood_function(const CSampleIDWeight &x);
	// returning value is the real log_likelihood calculated from target_model

protected:
	// Draw samples
	int EE_Draw(const CEESParameter &, CStorageHead &); 	// equi-energy draw
	int EE_Draw_RandomBlock(const CEESParameter &, CStorageHead &); 	// equi-energy draw using random blocks
public:
	double BurnIn(size_t burn_in_length);		// returns the maximum posteior during burn-in
	double BurnIn_RandomBlock(const CEESParameter &, CStorageHead &); 	// burn-in using random blocks
	bool Initialize(CStorageHead &, unsigned int start_bin, unsigned int end_bin, size_t desired_pool_size);	// Initialize model (setting values for current_sample) using bins indexed from start_bin through (including) end_bin. 
	bool InitializeWithBestSample(CStorageHead &storage, unsigned int start_bin, unsigned int end_bin); 		// Initialize model using the best sample in the bins indexed from start_bin to end_bin
	double Simulation_Within(const CEESParameter &, CStorageHead &storage, bool if_storage, const string &sample_file_name=string()); 	// Simulation within the same energy level (no jumping across levels). Returns the maximum posterior during simulation
	double Simulation_Within_RandomBlock(const CEESParameter &, CStorageHead &storage, bool if_storage, const string &sample_file_name=string());  
	double Simulation_Cross(const CEESParameter &, CStorageHead &storage, bool if_storage, const string &sample_file_name=string()); 	// Simulation across levels. Returns the maximum posterior during simulation 	
	double Simulation_Cross_RandomBlock(const CEESParameter &, CStorageHead &storage, bool if_storage, const string &sample_file_name=string()); 
	void HillClimb(size_t nSolution, CStorageHead &storage, const CEESParameter &parameter); 
	
///////////////////////////////////////////////////////////////////////////////////////////
// IO 
public:
	void SaveSampleToStorage(const CSampleIDWeight &sample, const CEESParameter &parameter, CStorageHead &storage); 
	bool InitializeFromFile(const string &file_name); 

///////////////////////////////////////////////////////////////////////////////////////////
// Construction & Destruction functions here
public:
	CEquiEnergyModel(); 
	CEquiEnergyModel(bool _if_bounded, unsigned int eL, double _h, double _t, const CSampleIDWeight & _x=CSampleIDWeight(), TStateModel *_model=NULL, CMetropolis *_metropolis =NULL, time_t _time=time(NULL)); 
	~CEquiEnergyModel(); 
	/* Functions
	*/
friend class MinusLogPosterior_NPSOL; 
};

#endif
