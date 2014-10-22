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
 public:
  TTimeSeries *target;	
  TStorage *storage;
  CEESParameter *parameter;

  int nParameters;
  int nNodes;
  int node;             

  int stage;
  TDraw current_sample; 
  TDenseVector lambda;
  TDenseVector LogIntegralKernel;
  TDenseVector ESS;
  double scale;

  // Metropolis jumping kernel
  TDenseMatrix OrthonormalDirections;
  TDenseVector SqrtDiagonal;

  // initial distributions
  TDenseMatrix InitialOrthonormalDirections;
  TDenseVector InitialSqrtDiagonal;
  TDenseVector InitialCenter;

  // diagnostics
  TDenseVector nEEJumps, nEEProposed, nMHJumps, nMHProposed, nSaved, nStarts, striation_count, striation_up, striation_down;

  // Constructors/Destructors
  CEquiEnergyModel(int _nNode, int _node, TTimeSeries *_target, CEESParameter *_parameter);
  virtual ~CEquiEnergyModel() { delete storage; };

  // Setup
  void WriteMHInfo(void);
  void ReadMHInfo(int Stage);
  void SetupComputeNode(int Stage);
  void SetupMasterNode(const TDenseVector &mean, const TDenseMatrix &variance);
  void SetupMasterNode(int Stage);

  void ComputeUnweightedMeanVariance(int stage, TDenseVector &mean, TDenseMatrix &variance);
  void ComputeWeightedMeanVariance(int stage, const TDenseVector &importance_weights, TDenseVector &mean, TDenseMatrix &variance);

  double Tune(bool use_EE, double mid=0.3, int period=10, int max_period=0, bool verbose=true);
  void SimulateEE(bool store, int G, int N, bool start_from_current_sample=false, bool verbose=true);
  void SimulateMH(bool store, int number_to_save, bool start_from_current_sample=false, bool verbose=true);

  void ImportanceWeightedDraw(TDraw &x);
  void EqualWeightedDraw(TDraw &x, int striation);
  double GetMetropolisProposal(TDraw &x);
  double LogDensity(TDraw &x, int stage);

  void ResetDiagnostics(int number_striations);
  void WriteDiagnosticComputeNode(void);
  void WriteDiagnosticMasterNode(int stage);

  double InitialDraw(TDenseVector &x);

  // seting up filename
  string MakeFilename(const string &id, int level=-1, int node=-1);
  void OpenFile(fstream &file, bool output_file, const string &id, int level=-1, int node=-1);

  // initialization
  void WriteParameters(void);
  void WriteInitialDistribution(void);
  void ReadInitialDistribution(void);
};

#endif
