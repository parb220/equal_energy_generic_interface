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

	// Draw samples
	int EE_Draw(size_t MH_thin); 	// equi-energy draw
public:
	double BurnIn(size_t burn_in_length);		// returns the maximum posteior during burn-in
	bool Initialize(size_t desired_pool_size, int level);	// Initialize model (setting values for current_sample) using bins indexed from start_bin through (including) end_bin. 
	bool InitializeWithBestSample(int level); 		// Initialize model using the best sample in the bins indexed from start_bin to end_bin
	bool InitializeWith_Kth_BestSample(size_t K, int level_index);
	bool Initialize_RandomlyPickFrom_K_BestSample(size_t K, int level_index); 
	bool k_means_clustering(size_t K, int level_index, vector<CSampleIDWeight> &centers) const; 
	bool InitializeFromFile(const string &file_name); 

	double Simulation_Within(bool if_storage, const string &sample_file_name=string()); 	// Simulation within the same energy level (no jumping across levels). Returns the maximum posterior during simulation
	double Simulation_Cross(bool if_storage, const string &sample_file_name=string()); 	// Simulation across levels. Returns the maximum posterior during simulation 	

//////////////////////////////////////////////////////////////////////////////////////////
// 	make a equi-energy jump
	bool MakeEquiEnergyJump(CSampleIDWeight &y_end, const CSampleIDWeight &y_initial) const; 

///////////////////////////////////////////////////////////////////////////////////////////// Random Block (not debugged yet)
	int EE_Draw_RandomBlock(size_t MH_thin); 	// equi-energy draw using random blocks
	double BurnIn_RandomBlock(size_t); 	// burn-in using random blocks
	double Simulation_Within_RandomBlock(bool if_storage, const string &sample_file_name=string());  
	double Simulation_Cross_RandomBlock(bool if_storage, const string &sample_file_name=string()); 

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
};

#endif
