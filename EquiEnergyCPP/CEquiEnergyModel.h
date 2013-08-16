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

class CMetropolis; 
class CEquiEnergyModel; 

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
 	unsigned int energy_level;		 
	double h_bound, t_bound;
//////////////////////////////////////////////////////////////
// Current sample holding the following 
// 	Sample: 
// 	ID:	unique id of the sample (by default is the difference
// 		between when it is drawn and when the program is started)
// 	Weight:	energy (log posterior)
	CSampleIDWeight current_sample; 
	CMetropolis *metropolis; 	// pointer to CMetropolis
	time_t timer_when_started; 
	
public:
	virtual double log_posterior_function(CSampleIDWeight &x)=0;
        virtual double log_likelihood_function(const CSampleIDWeight &x)=0;

	// Draw samples
	int EE_Draw(const CEESParameter &, CStorageHead &, size_t MH_thin); 	// equi-energy draw
	int EE_Draw_RandomBlock(const CEESParameter &, CStorageHead &, size_t MH_thin); 	// equi-energy draw using random blocks
public:
	double BurnIn(size_t burn_in_length);		// returns the maximum posteior during burn-in
	double BurnIn_RandomBlock(size_t); 	// burn-in using random blocks
	bool Initialize(CStorageHead &, unsigned int start_bin, unsigned int end_bin, size_t desired_pool_size);	// Initialize model (setting values for current_sample) using bins indexed from start_bin through (including) end_bin. 
	bool InitializeWithBestSample(CStorageHead &storage, unsigned int start_bin, unsigned int end_bin); 		// Initialize model using the best sample in the bins indexed from start_bin to end_bin
	bool InitializeWith_Kth_BestSample(unsigned int K, CStorageHead &storage, unsigned int start_bin, unsigned int end_bin);
	bool Initialize_RandomlyPickFrom_K_BestSample(size_t K, CStorageHead &storage, unsigned int start_bin, unsigned int end_bin); 
	bool InitializeFromFile(const string &file_name); 

	double Simulation_Within(const CEESParameter &, CStorageHead &storage, bool if_storage, const string &sample_file_name=string()); 	// Simulation within the same energy level (no jumping across levels). Returns the maximum posterior during simulation
	double Simulation_Within_RandomBlock(const CEESParameter &, CStorageHead &storage, bool if_storage, const string &sample_file_name=string());  
	double Simulation_Cross(const CEESParameter &, CStorageHead &storage, bool if_storage, const string &sample_file_name=string()); 	// Simulation across levels. Returns the maximum posterior during simulation 	
	double Simulation_Cross_RandomBlock(const CEESParameter &, CStorageHead &storage, bool if_storage, const string &sample_file_name=string()); 
	
///////////////////////////////////////////////////////////////////////////////////////////
// IO 
public:
	virtual void SaveSampleToStorage(const CSampleIDWeight &sample, unsigned int, CStorageHead &storage); 
	virtual void Take_Sample_Just_Drawn_From_Storage(const CSampleIDWeight &x); 

///////////////////////////////////////////////////////////////////////////////////////////
// Construction & Destruction functions here
public:
	CEquiEnergyModel(); 
	CEquiEnergyModel(bool _if_bounded, unsigned int eL, double _h, double _t, const CSampleIDWeight & _x=CSampleIDWeight(), CMetropolis *_metropolis =NULL, time_t _time=time(NULL)); 
	virtual ~CEquiEnergyModel() {} 
};

#endif
