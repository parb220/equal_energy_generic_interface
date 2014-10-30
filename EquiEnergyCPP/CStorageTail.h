#ifndef _CLASS_STORAGE_TAIL_
#define _CLASS_STORAGE_TAIL_

#include <vector>
#include <fstream>
#include "dw_dense_matrix.hpp"

using namespace std; 

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
   There are n_parameters + 4 elements in each row of the matrices in or out.
     (1) 0 : n_parameters-1 - parameters 
     (2) n_parameters       - log posterior
     (3) n_parameters+1     - log kernel 
     (4) n_parameters+2     - group index
     (5) n_parameters+3     - sequence index
*/
class TStorage
{
 protected:
  // common parameters
  int n_parameters;
  int size;
  string storage_directory;
  string run_id;
  int node;

  // output
  TDenseMatrix out;
  int idx_out;
  int stage_out;
  int file_number_out;

  // cumulative importance weights                        
  TDenseVector weights;     

  // striations
  int n_striations;
  vector<int> striation_offset;
  vector<int> striation_width;
  TDenseVector striation_boundary;

  // input
  TDenseMatrix in_weighted;
  int idx_in_weighted;
  vector<TDenseMatrix> in;
  vector<int> idx_in;

  void Reset(void);
  void WriteDraws(void);
  void FetchEqualWeightedDraws(int striation);
  void FetchImportanceWeightedDraws(void);

 public:  
  TStorage(int NumberParameters, const string &StorageDirectory, const string &RunId, int Node);
  ~TStorage() {};

  void Setup(int stage, int number_striations, double lambda);

  void FinalizeComputeNode(void);
  void FinalizeMasterNode(int stage);

  void ReadInfo(int stage, TDenseVector &log_posterior, TDenseVector &log_density);

  // output (storing draws)
  void AddDraw(const TDraw &x);

  // input (getting draws)
  void EqualWeightedDraw(TDraw &x, int striation);
  void ImportanceWeightedDraw(TDraw &x);

  // striations
  int NumberStriations(void) { return n_striations; };
  int Striation(double log_posterior);
  double UpperStriationBoundary(int striation);
  double LowerStriationBoundary(int striation);
  double WeightedStriationProbability(int striation);

  // creating files
  string MakeFullFileName(string id, int level=-1, int node=-1, int file_number=-1);
  void OpenFile(fstream &file, bool output_file, const string &id, int level=-1, int node=-1, int file_number=-1);

  void WriteASCIIDraws(int stage);

  void CheckClass(void);

  // Removing file
  void RemoveRecords(int stage); 
};



#endif
