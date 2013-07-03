#ifndef _CLASS_STATE_MODEL
#define _CLASS_STATE_MODEL

#include <vector>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <cstdio>
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CIndexBlockSizeScale.h"

extern "C" {
	#include "dw_switch.h"
}

using namespace std;

const int METROPOLIS_JUMP = 1; 
const int EQUI_ENERGY_JUMP = -1; 
const int NO_JUMP = 0; 

class CEquiEnergyModel
{
protected:
	time_t timer_when_started; 
public:
	TStateModel *target_model; 	// pointer to the TStateModel 
	CMetropolis *metropolis;	// pointer to metropolis

///////////////////////////////////////////////////////////////
// Parameters for equi-energy sampling
// 	energy_level: i
// 	energy_bound: h[i]
// 	temperature_bound: t[i]
	int energy_level;		 
	double h_bound, t_bound;

//////////////////////////////////////////////////////////////
// Current sample holding the following 
// 	Sample: 
// 	ID:	unique id of the sample (by default is the difference
// 		between when it is drawn and when the program is started)
// 	Weight:	energy (log posterior)
	CSampleIDWeight current_sample; 

///////////////////////////////////////////////////////////////
// Calculate log-posterior and log-likelihood of a given sample
// Depending on TStateModel, LogPosterior and LogLikelihood can be implemented 
// differently. 
public:
	double LogPosterior(const CSampleIDWeight &x)
	{
		return LogPosterior(x.data.vector, x.data.dim); 
	}
	double LogLikelihood(const CSampleIDWeight &x)
	{
		return LogLikelihood(x.data.vector, x.data.dim); 
	}
	double LogPosterior(const TDenseVector &x) 
	{
		return LogPosterior(x.vector, x.dim);
	}
	double LogLikelihood(const TDenseVector &x)
	{
		return LogLikelihood(x.vector, x.dim); 
	}
	virtual double LogPosterior(const double *x, int nx); 
	virtual double LogLikelihood(const double *x, int nx);

protected:
	// Draw samples
	int EE_Draw(const CEESParameter &, CStorageHead &); 	// equi-energy draw
	double BurnIn(size_t burn_in_length);		// returns the maximum posteior during burn-in
	bool Initialize(CStorageHead &, unsigned int start_bin, unsigned int end_bin, size_t desired_pool_size);	// Initialize model (setting values for current_sample) using bins indexed from start_bin through (including) end_bin. 
	bool InitializeWithBestSample(CStorageHead &storage, unsigned int start_bin, unsigned int end_bin); 		// Initialize model using the best sample in the bins indexed from start_bin to end_bin
	double Simulation_Within(const CEESParameter &, CStorageHead &storage, bool if_storage, FILE *, bool if_write_file); 	// Simulation within the same energy level (no jumping across levels). Returns the maximum posterior during simulation
	double Simulation_Cross(const CEESParameter &, CStorageHead &storage, bool if_storage, FILE *, bool if_write_file); 	// Simulation across levels. Returns the maximum posterior during simulation 
	
	
///////////////////////////////////////////////////////////////////////////////////////////
// IO 
public:
	void SaveSampleToStorage(const CSampleIDWeight &sample, const CEESParameter &parameter, CStorageHead &storage); 
	void SaveSampleToFile(const CSampleIDWeight &sample, const CEESParameter &parameter, FILE *file); 

protected:
	// Functions
	int EE_Draw_MultipleBlock(const double*, int *, int *, int, const CEESParameter &, CStorageHead &, const gsl_rng *); 	// equi-energy draw in multiple blocks
public:


	double BurnIn_MultipleBlock(const vector <vector<CIndexBlockSizeScale> > &, const CEESParameter &, int); 

	double Simulation_Within_MultipleBlock(const vector<vector<CIndexBlockSizeScale> > &, const CEESParameter &, CStorageHead &storage, bool if_storage, FILE *, bool if_write_file, const gsl_rng *);
	
	double Simulation_Cross_MultipleBlock(const vector<vector<CIndexBlockSizeScale> > &, const CEESParameter &, CStorageHead &storage, bool if_storage, FILE *, bool if_write_file, const gsl_rng *);


	bool LoadMetropolisParameter(string); 
///////////////////////////////////////////////////////////////////////////////////////////// Construction, deconstruction
public: 
	CEquiEnergyModel(TStateModel *_target=NULL, CMetropolis *_metropolis=NULL, int _level=0, double _h=0 double _t=1, const CSampleIDWeight &_x=CSampleIDWeight() ); 
	virtual ~CEquiEnergyModel();
	

	/* All Calibrate functions return a vector of thetaDim*2 containing the scales 
	// of each dimension during the two-passes training
	int Calibrate_Diag(vector<CIndexBlockSizeScale> &, vector<CIndexBlockSizeScale> &, double, double, int, int, int, int);
	int Calibrate_Variance_Initialize_Sample(vector<CIndexBlockSizeScale> &, vector<CIndexBlockSizeScale> &, string, double, double, int, int, int, int); 
	int Calibrate_Diag_Random_Block(vector<CIndexBlockSizeScale> &, const vector <CIndexBlockSizeScale> &, const vector<CIndexBlockSizeScale> &, double, int, int, const CEESParameter&); 
	int Calibrate_Variance_Initialize_Sample_Random_Block(vector<CIndexBlockSizeScale> &, const vector<CIndexBlockSizeScale> &, const vector<CIndexBlockSizeScale> &, string, double, int, int, const CEESParameter &); 
	*/

};

#endif
