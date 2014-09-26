#ifndef _CLASS_STATE_MODEL
#define _CLASS_STATE_MODEL

#include <vector>
#include <ctime>
#include <cstdio>
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "dw_time_series.hpp"

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

class CEquiEnergyModel
{
private:
	CEquiEnergyModel(const CEquiEnergyModel &);

public:
///////////////////////////////////////////////////////////////
	
public:
	TTimeSeries *target;	
	CStorageHead *storage; 
	CEESParameter *parameter;

	int nParameters;
	int nNodes;
 	int energy_level;
	CSampleIDWeight current_sample; 
	time_t timer_when_started; 
public:
	bool Initialize_WeightedSampling(int K, int level_index, vector<CSampleIDWeight> &starters) const; 

///////////////////////////////////////////////////////////////////////////////////////////
// IO 
public:
	virtual void SaveSampleToStorage(const CSampleIDWeight &sample); 
	virtual void Take_Sample_Just_Drawn_From_Storage(const CSampleIDWeight &x); 

///////////////////////////////////////////////////////////////////////////////////////////
// Construction & Destruction functions here
public:
	CEquiEnergyModel(int nNode, int _rank, TTimeSeries *_target, CEESParameter *_parameter);
	virtual ~CEquiEnergyModel() { delete storage; }; 

//////////////////////////////////////////////////////////////////////////////////////////
// Sampling and Tuning - added by DW 9/10/2014
public:
	TDenseVector K;
	TDenseVector LogKernelIntegral;
	double scale;
	double density_constant;
	int node;

	int nEEJumps, nEEProposed, nMHJumps, nMHProposed, nRestarts, nSaved, ring_changes_up, ring_changes_down;
	vector<int> ring_counts;

        double Tune(double mid=0.3, int period=10, int max_period=0, bool verbose=true);
	void SimulateEE(bool if_storage, const string &sample_file_name, int number_to_save, int reinitialize_factor, bool start_from_current_sample=false, bool verbose=true);
	void SimulateMH(bool store_internal, const string &external_filename, int number_to_save, bool verbose);
        double GetMetropolisProposal(CSampleIDWeight &x_new);

	// importance sampling
	vector<CSampleIDWeight> ImportanceSamples;
	int iImportanceSamples;
	bool GetImportanceSamples(void);
	bool ImportanceSamplePreviousLevel(CSampleIDWeight &x_new);

	double ESS;
	void WriteSimulationDiagnostic(int node);
	void WriteSimulationDiagnostic(void);

	/////////////////////////////////////////////////////////////////////////
	TDenseMatrix OrthonormalDirections;
	TDenseVector SqrtDiagonal;

	// initial distributions
	TDenseMatrix InitialOrthonormalDirections;
	TDenseVector InitialSqrtDiagonal;
	TDenseVector InitialCenter;

	double InitialDraw(TDenseVector &x);
	double LogInitialDensity(const TDenseVector &x);
	double AnalyzeInitialDraws(void);

	// seting up filename
	string MakeFilename(const string &id, int level=-1, int node=-1);
	void OpenFile(fstream &file, bool output_file, const string &id, int level=-1, int node=-1);

	// initialization
	void SetupFromPreviousLevel(int level);
	void WriteInitializationFile(void);
	void ReadInitializationFile(int level);
	void WriteScale(int level, int node, double s);
	void WriteScale(int level, double s);
	double ConsolidateScales(int level);
	void WriteParameters(void);

};

double EffectiveSampleSize(TDenseVector log_kernel_target, TDenseVector log_density_proposal);
double ComputeLogKernelIntegral(TDenseVector log_kernel_target, TDenseVector log_density_proposal);

#endif
