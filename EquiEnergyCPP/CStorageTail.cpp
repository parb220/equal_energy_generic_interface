
#include <sys/types.h> 
#include <sys/stat.h>
#include <sys/file.h>
#include <errno.h>
#include <fstream>
#include <iomanip>

#include "CStorageTail.h"
#include "dw_rand.h"
#include "dw_math.h"

using namespace std;

#define N_DRAWS 1024

/*
  Assume that v[i-1] <= v[i] for 0 < i < v.dim.  If
     
  u < v[0]             returns 0
  v[i-1] <= u < v[i]   returns i
  v[v.dim-1] <= u      returns v.dim
*/
static int StrictLeastUpperBound(const TDenseVector &v, double u)
{
  // v[lo] <= u < v[hi]
  int lo=0,  hi=v.dim-1, mid;
  if (u < v.vector[lo]) return lo;
  if (v.vector[hi] <= u) return hi+1;
  while (hi - lo > 7)
    {
      mid=(hi+lo)/2;
      if (u < v.vector[mid])
  	hi=mid;
      else
  	lo=mid;
    }
  if (u < v.vector[++lo]) return lo;
  if (u < v.vector[++lo]) return lo;
  if (u < v.vector[++lo]) return lo;
  if (u < v.vector[++lo]) return lo;
  if (u < v.vector[++lo]) return lo;
  if (u < v.vector[++lo]) return lo;
  return hi;
}

static TDenseVector ComputeCumulativeImportanceWeights(double lambda, const TDenseVector &log_posterior, const TDenseVector &log_density)
{
  double sum=log_posterior.vector[0]*lambda - log_density.vector[0];
  for (int i=1; i < log_posterior.dim; i++)
    sum=AddLogs(sum,log_posterior.vector[i]*lambda - log_density.vector[i]);

  TDenseVector cumulative_importance_weights(log_posterior.dim);
  cumulative_importance_weights.vector[0]=exp(log_posterior.vector[0]*lambda - log_density.vector[0] - sum);
  for (int i=1; i < log_posterior.dim; i++)
    cumulative_importance_weights.vector[i]=cumulative_importance_weights.vector[i-1] + exp(log_posterior.vector[i]*lambda - log_density.vector[i] - sum);
 
  return cumulative_importance_weights;
}

/////////////////////////////////////////////////////////////////////////////////
/// Begin TStorage
/////////////////////////////////////////////////////////////////////////////////
TStorage::TStorage(int NumberParameters, const string &StorageDirectory, const string &RunId, int Node) 
  : n_parameters(NumberParameters), size(NumberParameters+4), storage_directory(StorageDirectory), run_id(RunId), node(Node), 
    out(8192,size,false), idx_out(0), stage_out(-1), file_number_out(0), weights(0),
    n_striations(0), striation_offset(0), striation_width(0), striation_boundary(0), 
    in_weighted(0,size,false), idx_in_weighted(0), in(0), idx_in(0)
{
  // ensure there is a slash at the end of storage directory
  if (!storage_directory.empty() && (storage_directory.at(storage_directory.length()-1) != '/')) storage_directory.push_back('/');
  storage_directory=storage_directory + run_id;
  storage_directory.push_back('/');

  // if master node, make subdirectory for run
  if (node == 0)
    {
      errno = 0; 
      int status = mkdir(storage_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (status !=0 && errno != EEXIST)
	{
	  cerr << "Unable to create directory " << storage_directory << endl;
	  abort();
	}
    }
}

void TStorage::Reset(void)
{
  // output
  idx_out=0;
  stage_out=-1;
  file_number_out=0;

  // cumulative importance weights
  weights.Resize(0);

  // striations
  n_striations=0;
  striation_offset.resize(0);
  striation_width.resize(0);
  striation_boundary.Resize(0);

  // input 
  in_weighted.Resize(0,size);
  idx_in_weighted=0;
  in.resize(0);
  idx_in.resize(0);
}

/*
   Writes any unsaved draws to disk and resets all parameters to default values.
*/
void TStorage::FinalizeComputeNode(void)
{
  // write any unsaved draws to disk
  WriteDraws();

  // set all values to default values
  Reset();
}

/*
  Consolidates all draws into a single file, which is sorted in ascending order 
  by the log posterior.  The original files are
  deleted.

  Creates info file consisting of the stage, the vector of log posteriors,
  and the vector of log densities.

  Currently all draws are loaded into memory, sorted, and then written back to 
  disk.
*/
void TStorage::FinalizeMasterNode(int stage)
{
  string filename;
  fstream file;
  TDenseMatrix draws(0,size,false), matrix;
  int node=1;
  while (1)
    {
      int file_number=-1;
      while (1)
	{
	  filename=MakeFullFileName("draws",stage,node,file_number+1);
	  file.open(filename.c_str(), ios::in|ios::binary);
	  if (!file.is_open()) break;
	  matrix.ReadBinary(file);
	  draws=VCat(draws,matrix);
	  file.close();
	  remove(filename.c_str()); 
	  file_number++;
	}
      if (file_number == -1) break;
      node++;
    }

  if (node == 1)
    {
      cerr << "No draws files to consolidate" << endl;
      abort();
    }

  // write draws file
  draws.SortRows(n_parameters);
  OpenFile(file,true,"draws",stage);
  draws.WriteBinary(file);
  file.close();

  // write info file
  TDenseVector log_posterior=draws.ColumnVector(n_parameters), log_density=draws.ColumnVector(n_parameters+1);
  OpenFile(file,true,"info",stage);
  file.write((char*)(&stage),sizeof(int));
  log_posterior.WriteBinary(file);
  log_density.WriteBinary(file);
  file.close();
}

void TStorage::RemoveRecords(int stage)
{
	string filename=MakeFullFileName("draws",stage);
	remove(filename.c_str()); 
	filename = MakeFullFileName("info", stage); 
	remove(filename.c_str()); 
	filename = MakeFullFileName("Diagnostic", stage); 
	remove(filename.c_str());
	filename = MakeFullFileName("MH-info", stage); 
	remove(filename.c_str()); 
}

void TStorage::ReadInfo(int stage, TDenseVector &log_posterior, TDenseVector &log_density)
{
  fstream file;
  OpenFile(file, false, "info", stage);
  int info_stage;
  file.read((char*)(&info_stage),sizeof(int));
  if (info_stage != stage)
    {
      cerr << "TStorage::ReadInfo(): name and internal stored stage do not agree" << endl;
      abort();
    }
  log_posterior.ReadBinary(file);
  log_density.ReadBinary(file);
  file.close();
}

void TStorage::Setup(int stage, int number_striations, double lambda)
{
  // reset all fields
  Reset();

  // output
  idx_out=0;
  stage_out=stage;
  file_number_out=0;

  if (stage > 0)
    {
      // cumulative importance weights
      TDenseVector log_posterior, log_density;
      ReadInfo(stage-1,log_posterior,log_density);
      weights=ComputeCumulativeImportanceWeights(lambda,log_posterior,log_density);

      // striations
      n_striations=number_striations;
      striation_offset.resize(n_striations);
      striation_width.resize(n_striations);
      striation_boundary.Resize(n_striations-1);

      // striation offset and width
      for (int m=1+(log_posterior.dim)/number_striations, hi=log_posterior.dim, i=n_striations-1; i >= 0; i--)
	{
	  if (hi == 0)
	    {
	      cerr << "SetupForInput(): No draws in some striations" << endl;
	      abort();
	    }
	  int j = (hi > m) ? hi-m : 0;
	  for ( ; (j > 0) && (log_posterior.vector[j] == log_posterior.vector[j-1]); j--);
	  striation_width[i]=hi-j;
	  striation_offset[i]=hi=j;
	}

      // striation boundary
      for (int i=0; i < striation_boundary.dim; i++)
	striation_boundary[i]=log_posterior.vector[striation_offset[i+1]];

      // input
      in.resize(n_striations);
      idx_in.resize(n_striations);
      for (int i=n_striations-1; i >= 0; i--) idx_in[i]=0;
      idx_in_weighted=0;
    }
}

void TStorage::WriteDraws(void)
{
  if (idx_out > 0)
    {
      fstream file;
      OpenFile(file,true,"draws",stage_out,node,file_number_out++);
      out.WriteBinary(file,0,idx_out-1,0,size-1);
      file.close();
      idx_out=0;
    }
}

void TStorage::FetchEqualWeightedDraws(int striation)
{
  // open file
  fstream file_in;
  OpenFile(file_in,false,"draws",stage_out-1);

  in[striation].UniqueMemory(N_DRAWS,size,false);
  for (int i=0; i < N_DRAWS; i++)
    {
      // get index
      int idx=floor(striation_width[striation]*dw_uniform_rnd());
      if (idx >= striation_width[striation]) idx=striation_width[striation]-1;
      idx+=striation_offset[striation];

      // seek and read
      file_in.seekg((2*sizeof(int)+sizeof(bool))+idx*size*sizeof(double),ios::beg);
      file_in.read((char*)(in[striation].matrix+i*size),size*sizeof(double));
    }
  file_in.close();
  idx_in[striation]=N_DRAWS;
}

void TStorage::FetchImportanceWeightedDraws(void)
{
  // open file
  fstream file_in;
  OpenFile(file_in,false,"draws",stage_out-1);

  in_weighted.UniqueMemory(N_DRAWS,size,false);
  for (int i=0; i < N_DRAWS; i++)
    {
      // get index
      int idx=StrictLeastUpperBound(weights,dw_uniform_rnd()) - 1;

      // seek and read
      file_in.seekg((2*sizeof(int)+sizeof(bool))+idx*size*sizeof(double),ios::beg);
      file_in.read((char*)(in_weighted.matrix+i*size),size*sizeof(double));
    }
  file_in.close();
  idx_in_weighted=N_DRAWS;
}

void TStorage::EqualWeightedDraw(TDraw &x, int striation)
{
  if ((striation < 0) || (striation >= n_striations))
    {
      cerr << "TStorage::EqualWeightedDraw(): invalid striation index" << endl;
      abort();
    }
  if (idx_in[striation] == 0) FetchEqualWeightedDraws(striation);
  int idx=--idx_in[striation];
  x.Set(in[striation].RowVector(idx,0,n_parameters-1), in[striation](idx,n_parameters),in[striation](idx,n_parameters+1), 
	in[striation](idx,n_parameters+2), in[striation](idx,n_parameters+3));
}

void TStorage::ImportanceWeightedDraw(TDraw &x)
{
  if (idx_in_weighted == 0) FetchImportanceWeightedDraws();
  idx_in_weighted--;
  x.Set(in_weighted.RowVector(idx_in_weighted,0,n_parameters-1), in_weighted(idx_in_weighted,n_parameters), in_weighted(idx_in_weighted,n_parameters+1), 
	in_weighted(idx_in_weighted,n_parameters+2), in_weighted(idx_in_weighted,n_parameters+3));
}

void TStorage::AddDraw(const TDraw &x)
{
  if (idx_out == out.rows) WriteDraws();
  out.InsertRowMatrix(idx_out,0,x.Parameters());
  out(idx_out,n_parameters)=x.LogPosterior();
  out(idx_out,n_parameters+1)=x.LogKernel();
  out(idx_out,n_parameters+2)=x.Group();
  out(idx_out,n_parameters+3)=x.Index();
  idx_out++;
}

string TStorage::MakeFullFileName(string id, int level, int node, int file_number)
{
  stringstream filename;
  filename << storage_directory << run_id << '.' << id;
  if (level >= 0)
    {
      filename << '.' << level;
      if (node >= 0)
	{
	  filename << '.' << node;
	  if (file_number >= 0)
	    filename  << '.' << file_number;
	}
    }
  return filename.str();
}

void TStorage::OpenFile(fstream &file, bool output_file, const string &id, int level, int node, int file_number)
{
  string filename=MakeFullFileName(id,level,node,file_number);
  file.open(filename.c_str(), output_file ? ios::out|ios::binary : ios::in|ios::binary );
  if (!file.is_open())
    {
      cerr << "TStorage::OpenFile(): Unable to open " << filename << endl;
      abort();
    }
}

int TStorage::Striation(double log_posterior) 
{
  return StrictLeastUpperBound(striation_boundary, log_posterior);
}

double TStorage::UpperStriationBoundary(int striation) 
{ 
  if ((striation < 0) || (striation >= n_striations))
    {
      cerr << "TStorage::UpperStriationBoundary(): invalid striation index" << endl;
      abort();
    }
  return (striation < (int)striation_boundary.dim) ? striation_boundary[striation] : 1E300;
}

double TStorage::LowerStriationBoundary(int striation) 
{ 
  if ((striation < 0) || (striation >= n_striations))
    {
      cerr << "TStorage::LowerStriationBoundary(): invalid striation index" << endl;
      abort();
    }
  return (striation == 0) ? -1.0E300 : striation_boundary[striation-1];
}

double TStorage::WeightedStriationProbability(int striation) 
{ 
  if ((striation < 0) || (striation >= n_striations))
    {
      cerr << "TStorage::WeightedStriationProbability(): invalid striation index" << endl;
      abort();
    }
  return weights[striation_offset[striation]+striation_width[striation]-1]-weights[striation_offset[striation]];
}
/////////////////////////////////////////////////////////////////////////////////
/// End TStorage
/////////////////////////////////////////////////////////////////////////////////

void TStorage::CheckClass(void)
{
  cout << "TStorage::CheckClass() - stage" << stage_out << endl;

  if (stage_out == 0) return;

  if (NumberStriations() <= 0)
    {
      cerr << "class not setup for input" << endl;
      abort();
    }

  // check weights
  if (weights(0) <= 0)
    {
      cerr << "non-positive weight: " << weights(0) << endl;
      abort();
    }
  for (int i=1; i < weights.dim; i++)
    if (weights(i) <= weights(i-1))
      {
	cerr << "cumulative weights not ascending" << weights(i) << " " << weights(i+1) << endl;
	abort();
      }

  // check weighted draws
  cout << "Checking weighted draws" << endl;
  int ndraws=100000;
  TDraw x;
  TDenseVector c(NumberStriations(),0.0);
  for (int i=0; i < ndraws; i++)
    {
      ImportanceWeightedDraw(x);
      c(Striation(x.LogPosterior()))+=1;
    }
  for (int i=0; i < NumberStriations(); i++)
    cout << c(i)/(double)ndraws << "(" << WeightedStriationProbability(i) << ")" << endl;
  cout << endl;

  // check striation assignment
  cout << "Checking striation assignment" << endl;
  fstream file;
  int rows, s;
  TDenseVector z(size);
  OpenFile(file,false,"draws",stage_out-1);
  file.read((char*)(&rows),sizeof(int));
  file.seekg(2*sizeof(int)+sizeof(bool),ios::beg);
  for (int i=0; i < rows; i++)
    {
      file.read((char*)(z.vector),size*sizeof(double));
      s=Striation(z(n_parameters));
      if ((z(n_parameters) < LowerStriationBoundary(s)) || (z(n_parameters) >= UpperStriationBoundary(s)))
	cerr << "Invalid striation assignment " << LowerStriationBoundary(s) << " <= " << z(n_parameters) << " < " << UpperStriationBoundary(s) << endl;
    }
  cout << "Done checking striation assignment" << endl;

  // check unweighted draws
  for (int i=0; i < NumberStriations(); i++)
    {
      cout << "Checking unweighted draws - striation " << i << endl;
      for (int k=0; k < ndraws; k++)
	{
	  EqualWeightedDraw(x,i);
	  if (Striation(x.LogPosterior()) != i)
	    cerr << "draws not from striation" << endl;
	}
    }

  // speed tests
  cout << "Speed tests" << endl;

  time_t start_time=time((time_t*)NULL);
  for (int i=0; i < ndraws; i++)
    ImportanceWeightedDraw(x);
  cout << "Ellapsed time - weighted draws: " << difftime(time((time_t*)NULL),start_time) << " seconds " << endl;

  start_time=time((time_t*)NULL);
  for (int i=0; i < ndraws; i++)
    EqualWeightedDraw(x,0);
  cout << "Ellapsed time - unweighted draws(0): " << difftime(time((time_t*)NULL),start_time) << " seconds " << endl;

  start_time=time((time_t*)NULL);
  for (int i=0; i < ndraws; i++)
    EqualWeightedDraw(x,NumberStriations()/2);
  cout << "Ellapsed time - unweighted draws(" << NumberStriations()/2 << "): " << difftime(time((time_t*)NULL),start_time) << " seconds " << endl;

  start_time=time((time_t*)NULL);
  for (int i=0; i < ndraws; i++)
    EqualWeightedDraw(x,NumberStriations()-1);
  cout << "Ellapsed time - unweighted draws(" << NumberStriations()-1 << "): " << difftime(time((time_t*)NULL),start_time) << " seconds " << endl;
}
