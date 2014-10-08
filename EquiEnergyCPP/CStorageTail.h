#ifndef _CLASS_STORAGE_TAIL
#define _CLASS_STORAGE_TAIL

#include <vector>
#include <fstream>
#include "dw_dense_matrix.hpp"
#include "CSampleIDWeight.h"

int UpperBound(const TDenseVector &v, double u);

class TDraw
{
 protected:
  double log_posterior;
  double log_kernel;
  int group;
  int index;
  TDenseVector parameters;

 public:
  TDraw(void) {};
  ~TDraw() {};

  void Set(const TDenseVector &x, double Log_Posterior, double Log_Kernel, int Group, int Index) { parameters=x; log_posterior=Log_Posterior; log_kernel=Log_Kernel; group=Group; index=Index; };
  void Set(const TDenseVector &x, double Log_Posterior, double Log_Kernel) { parameters=x; log_posterior=Log_Posterior; log_kernel=Log_Kernel; group=-1; index=-1; };
  void SetInfo(double Group, double Index) { group=Group; index=Index; };
  void SetLogKernel(double Log_Kernel) { log_kernel=Log_Kernel; };

  const TDenseVector& Parameters(void) const { return parameters; };
  double LogPosterior(void) const { return log_posterior; };
  double LogKernel(void) const { return log_kernel; };
  int Group(void) const { return group; };
  int Index(void) const { return index; };
};


/*
   There are n_parameters + 2 elements in each row of the matrices in or out.
     (1) 0 : n_parameters-1 - parameters 
     (2) n_parameters       - log posterior
     (3) n_parameters+1     - log kernel 
     (4) n_parameters+2     - group index
     (5) n_parameters+3     - sequence index
*/
class TStorage
{
  //protected:
 public:
  int n_parameters;
  int size;
  string storage_directory;
  string run_id;

  TDenseMatrix out;
  int idx_out;
  int node_out;
  int file_number_out;
  int level_out;

  TDenseMatrix in;
  TDenseVector cumulative_importance_weights;
  TDenseVector importance_weights;
  vector<int> ring_offset;
  vector<int> ring_width;
  TDenseVector ring_boundary;
  double K_level_in_minus_one;
  double LogIntegralKernel;
  int level_in;

  void SetupForInput(int level, int number_rings);
  void WriteDraws(void);
  void SetupWeights(void);
  void SetTemperature(double K_max, double min_ess);

 public:
  TStorage(int NumberParameters, const string &StorageDirectory, const string &RunId);
  ~TStorage() {};

  // output (storing draws)
  void SetupForOutput(int level, int node);
  void AddDraw(const TDraw &x);
  void FlushOutput(void) { WriteDraws(); level_out=-1; };

  // input (getting draws)
  void SetupForInput(int level, int number_rings, double temperature);
  void EqualWeightedDraw(TDraw &x, int ring);
  void ImportanceWeightedDraw(TDraw &x);
  double SetupForInput(int level, int number_rings, double K_max, double min_ess);
  void FreeInput(void);

  void ConsolidateDraws(int level_in, int level_out);
  void DeleteDraws(int level);

  int NumberRings(void) { return ring_offset.size(); };
  int Ring(double posterior) { return UpperBound(ring_boundary,posterior); };

  double ComputeEffectiveSampleSize(void);
  double ComputeLogIntegralKernel(void) { return LogIntegralKernel; };
  void ComputeWeightedVarianceDecomposition(TDenseMatrix &U, TDenseVector &d);
  void ComputeVarianceDecomposition(TDenseMatrix &U, TDenseVector &d);
  void ComputeMean(TDenseVector &m);
  double ComputeWeightedStriationProbability(int ring) { return cumulative_importance_weights[ring_offset[ring]+ring_width[ring]-1]-cumulative_importance_weights[ring_offset[ring]]; };

  string MakeFullFileName(string id, int level=-1, int node=-1, int file_number=-1);
  string OpenFile(fstream &file, bool output_file, const string &id, int level=-1, int node=-1, int file_number=-1);

  double UpperRingBoundary(int ring) { return (ring < (int)ring_boundary.dim) ? ring_boundary[ring] : 1E300; };
  double LowerRingBoundary(int ring) { return (ring == 0) ? -1.0E300 : ring_boundary[ring-1]; };

  void PrintWeights(void);
};



#endif
