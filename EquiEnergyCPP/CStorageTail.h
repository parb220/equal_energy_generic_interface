#ifndef _CLASS_STORAGE_TAIL_
#define _CLASS_STORAGE_TAIL_

#include <vector>
#include <fstream>
#include "dw_dense_matrix.hpp"

#define LOG_WEIGHTS(log_prior, log_likelihood, log_density, lambda) ((log_prior + log_likelihood)*lambda - log_density)
#define LOG_DENSITY(log_prior, log_likelihood, log_integral, lambda) ((log_prior + log_likelihood)*lambda - log_integral)

//#define LOG_WEIGHTS(log_prior, log_likelihood, log_density, lambda) (log_prior + log_likelihood*lambda - log_density)
//#define LOG_DENSITY(log_prior, log_likelihood, log_integral, lambda) (log_prior + log_likelihood*lambda - log_integral)

using namespace std; 

class TDraw
{
 protected:
  double log_prior;
  double log_likelihood;
  double log_density;
  double log_density_previous;
  int group;
  int index;
  TDenseVector parameters;

 public:
  TDraw(void) {};
  ~TDraw() {};

  TDraw& operator=(const TDraw &x) 
    { log_prior=x.log_prior, log_likelihood=x.log_likelihood; log_density=x.log_density; log_density_previous=x.log_density_previous; group=x.group; index=x.index; parameters=x.parameters; return *this; };

  void Set(const TDenseVector &x, double Log_Prior, double Log_Likelihood, double Log_Density, int Group, int Index) 
    { parameters=x; log_prior=Log_Prior; log_likelihood=Log_Likelihood; log_density=Log_Density; log_density_previous=0.0; group=Group; index=Index; };
  void Set(const TDenseVector &x, double Log_Prior, double Log_Likelihood) 
    { parameters=x; log_prior=Log_Prior; log_likelihood=Log_Likelihood; log_density=log_density_previous=0.0; group=index=-1; };
  void SetInfo(double Group, double Index) { group=Group; index=Index; };
  void SetLogDensity(double Log_Density) { log_density=Log_Density; };
  void SetLogDensityPrevious(double Log_Density_Previous) { log_density_previous=Log_Density_Previous; };

  const TDenseVector& Parameters(void) const { return parameters; };
  double LogPrior(void) const { return log_prior; };
  double LogLikelihood(void) const { return log_likelihood; };
  double LogPosterior(void) const { return log_prior+log_likelihood; };
  double LogDensity(void) const { return log_density; };
  double LogDensityPrevious(void) const { return log_density_previous; };
  int Group(void) const { return group; };
  int Index(void) const { return index; };
};


/*
   There are n_parameters + 5 elements in each row of the matrices in or out.
     (1) 0 : n_parameters-1 - parameters 
     (2) n_parameters       - log prior
     (3) n_parameters+1     - log likelihood
     (4) n_parameters+2     - log kernel 
     (5) n_parameters+3     - group index
     (6) n_parameters+4     - sequence index
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

  void CumulativeImportanceWeights(double lambda, const TDenseVector &log_prior, const TDenseVector &log_likelihood, const TDenseVector &log_density);
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

  void ReadInfo(int stage, TDenseVector &log_prior, TDenseVector &log_likelihood, TDenseVector &log_density);

  // output (storing draws)
  void AddDraw(const TDraw &x);

  // input (getting draws)
  void EqualWeightedDraw(TDraw &x, int striation);
  void ImportanceWeightedDraw(TDraw &x);

  // striations
  int NumberStriations(void) const { return n_striations; };
  int Striation(const TDraw &x);
  double UpperStriationBoundary(int striation);
  double LowerStriationBoundary(int striation);
  double WeightedStriationProbability(int striation);

  // creating files
  string MakeFullFileName(string id, int level=-1, int node=-1, int file_number=-1);
  void OpenFile(fstream &file, bool output_file, const string &id, int level=-1, int node=-1, int file_number=-1);

  void WriteASCII_draws(int stage);
  void WriteASCII_info(int stage);

  void ErrorClass(void);
  void CheckClass(void);
};



#endif
