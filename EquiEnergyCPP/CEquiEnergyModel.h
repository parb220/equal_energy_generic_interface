#ifndef _CLASS_STATE_MODEL
#define _CLASS_STATE_MODEL

#include <vector>
#include <ctime>
#include <cstdio>
#include "CEESParameter.h"
#include "CStorageTail.h"
#include "dw_time_series.hpp"

using namespace std;

class CEquiEnergyModel
{
private:
	CEquiEnergyModel(const CEquiEnergyModel &);

public:
///////////////////////////////////////////////////////////////
	
public:
	TTimeSeries *target;	
	//CStorageHead *storage; 
	TStorage *storage;
	CEESParameter *parameter;

	int nParameters;
	int nNodes;
	int node;               // used by SimulateEE() to set group index
 	int energy_level;
	TDraw current_sample; 

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

	//double density_constant;

	int nMHJumps, nMHProposed, nRestarts, nSaved;
	vector<int> nEEJumps, nEEProposed;
	vector<int> ring_counts, ring_changes_up, ring_changes_down;

	void Setup(int level, double Kplus);
        double Tune(bool use_EE, double mid=0.3, int period=10, int max_period=0, bool verbose=true);
	void SimulateEE(bool store, int G, int N, bool start_from_current_sample=false, bool verbose=true);
	void SimulateMH(bool store, int number_to_save, bool start_from_current_sample=false, bool verbose=true);
        double GetMetropolisProposal(TDraw &x_new);

	TDenseVector ESS;
	void WriteSimulationDiagnostic(int level, int node);
	void WriteSimulationDiagnostic(int level);

	/////////////////////////////////////////////////////////////////////////
	TDenseMatrix OrthonormalDirections;
	TDenseVector SqrtDiagonal;

	// initial distributions
	TDenseMatrix InitialOrthonormalDirections;
	TDenseVector InitialSqrtDiagonal;
	TDenseVector InitialCenter;

	double InitialDraw(TDenseVector &x);
	double LogInitialDensity(const TDenseVector &x);

	// seting up filename
	string MakeFilename(const string &id, int level=-1, int node=-1);
	void OpenFile(fstream &file, bool output_file, const string &id, int level=-1, int node=-1);

	// initialization
	void WriteParameters(void);
	void WriteInitialDistribution(void);
	void ReadInitialDistribution(void);
};

#endif
