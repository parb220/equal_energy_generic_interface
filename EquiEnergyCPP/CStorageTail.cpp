
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

/*
   Assume that v[i-1] <= v[i] for 0 < i < v.dim.  If
     
                   u < v[0]   returns 0
         v[i-1] <= u < v[i]   returns i
     v[v.dim-1] <= u          returns v.dim
*/
int UpperBound(const TDenseVector &v, double u)
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


/////////////////////////////////////////////////////////////////////////////////
/// Begin TStorage
/////////////////////////////////////////////////////////////////////////////////
TStorage::TStorage(int NumberParameters, const string &StorageDirectory, const string &RunId) 
  : n_parameters(NumberParameters), size(NumberParameters+4), storage_directory(StorageDirectory), 
    run_id(RunId), out(8192,size,false), level_out(-1), level_in(-1) 
{
  // ensure there is a slash at the end of storage directory
  if (!storage_directory.empty() && (storage_directory.at(storage_directory.length()-1) != '/'))
    storage_directory.push_back('/');

  // make subdirectory for run
  storage_directory=storage_directory + run_id;
  errno = 0; 
  int status = mkdir(storage_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (status !=0 && errno != EEXIST)
    {
      cerr << "Unable to create directory " << storage_directory << endl;
      abort();
    }
  storage_directory.push_back('/');
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

string TStorage::OpenFile(fstream &file, bool output_file, const string &id, int level, int node, int file_number)
{
  string filename=MakeFullFileName(id,level,node,file_number);
  file.open(filename.c_str(), output_file ? ios::out : ios::in );
  if (!file.is_open())
    {
      cerr << "Unable to open " << filename << endl;
      abort();
    }
  return filename;
}

/*
  Sets up importance_weights and cumulative_importance_weights.
*/
void TStorage::SetupWeights(void)
{
  double sum=in(0,n_parameters)/K_level_in_minus_one - in(0,n_parameters+1);
  for (int i=1; i < in.rows; i++)
    sum=AddLogs(sum,in(i,n_parameters)/K_level_in_minus_one - in(i,n_parameters+1));

  importance_weights.Resize(in.rows);
  cumulative_importance_weights.Resize(in.rows);
  cumulative_importance_weights.vector[0]=importance_weights.vector[0]=exp(in(0,n_parameters)/K_level_in_minus_one - in(0,n_parameters+1) - sum);
  for (int i=1; i < in.rows; i++)
    {
      importance_weights.vector[i]=exp(in(i,n_parameters)/K_level_in_minus_one - in(i,n_parameters+1) - sum);
      cumulative_importance_weights.vector[i]=cumulative_importance_weights.vector[i-1] + importance_weights.vector[i];
    }

  LogIntegralKernel=sum - log(in.rows);
}

static double EffectiveSampleSize(TDenseVector log_kernel_target, TDenseVector log_density_proposal)
{
  double sum=log_kernel_target.vector[0] - log_density_proposal.vector[0];
  for (int i=1; i < log_kernel_target.dim; i++)
    sum=AddLogs(sum,log_kernel_target.vector[i] - log_density_proposal.vector[i]);

  double sum2=0;
  for (int i=0; i < log_kernel_target.dim; i++)
    {
      double weight=exp(log_kernel_target.vector[i] - log_density_proposal.vector[i] - sum);
      sum2+=weight*weight;
    }

  return 1/sum2;
}

double TStorage::ComputeEffectiveSampleSize(void)
{
  double inv_ess=0.0;
  for (int i=importance_weights.dim-1; i >= 0; i--) 
    inv_ess+=importance_weights.vector[i] * importance_weights.vector[i];
  return 1.0/inv_ess;
}

// double LogKernelIntegral(TDenseVector log_kernel_target, TDenseVector log_density_proposal)
// {
//   double sum=log_kernel_target.vector[0] - log_density_proposal.vector[0];
//   for (int i=log_kernel_target.dim-1; i > 0; i--)
//     sum=AddLogs(sum,log_kernel_target.vector[i] - log_density_proposal.vector[i]);
//   return sum - log(log_kernel_target.dim);
// }

// double TStorage::ComputeLogIntegralKernel(void)
// {
//   TDenseVector log_kernel_target=(1.0/K_level_in_minus_one)*in.ColumnVector(n_parameters);
//   TDenseVector log_kernel_proposal=in.ColumnVector(n_parameters+1);
//   return LogKernelIntegral(log_kernel_target,log_kernel_proposal);
// }

void TStorage::SetTemperature(double K_max, double min_ess)
{
  if (level_in > 1)
    {
      double ess, K_min=(K_max/2 <= 1.0) ? 1.0 : K_max/2;
      TDenseVector log_posterior=in.ColumnVector(n_parameters);
      TDenseVector log_kernel_proposal=in.ColumnVector(n_parameters+1);
      TDenseVector log_kernel_target;

      // bracket min_ess
      do
	{
	  log_kernel_target=(1.0/K_min)*log_posterior;
	  ess=EffectiveSampleSize(log_kernel_target,log_kernel_proposal);
	  if (ess >= min_ess)  K_min=(K_min/2 <= 1.0) ? 1.0 : K_min/2;
	} 
      while ((K_min > 1) && (ess >= min_ess));

      // bracketed?
      if (ess < min_ess)
	{
	  for (int i=0; i < 10; i++)
	    {
	      K_level_in_minus_one=0.5*(K_max + K_min);
	      log_kernel_target=(1.0/K_level_in_minus_one)*log_posterior;
	      ess=EffectiveSampleSize(log_kernel_target,log_kernel_proposal);

	      if (ess >= min_ess) 
		K_max=K_level_in_minus_one;
	      else
		K_min=K_level_in_minus_one;
	    }
	  K_level_in_minus_one=0.5*(K_max + K_min);
	}
      else
	K_level_in_minus_one=1.0;
    }
  else
    K_level_in_minus_one=1.0;
}

void TStorage::ComputeWeightedVarianceDecomposition(TDenseMatrix &U, TDenseVector &d)
{
  // Compute variance matrix
  TDenseVector v, m(n_parameters,0.0);
  TDenseMatrix variance(n_parameters,n_parameters,0.0);
  for (int i=in.rows-1; i >= 0; i--)
    {
      v=in.RowVector(i,0,n_parameters-1);
      m+=importance_weights.vector[i]*v;
      variance+=OuterProduct(importance_weights.vector[i]*v,v);
    }
  variance=variance - OuterProduct(m,m);

  // Compute decomposition
  Eig(d,U,0.5*(variance + Transpose(variance)));
  for (int j=d.dim-1; j >= 0; j--)  d.vector[j]=sqrt(d.vector[j]);
}

void TStorage::ComputeVarianceDecomposition(TDenseMatrix &U, TDenseVector &d)
{
  // Compute variance matrix
  TDenseVector v, m(n_parameters,0.0);
  TDenseMatrix variance(n_parameters,n_parameters,0.0);
  for (int i=in.rows-1; i >= 0; i--)
    {
      v=in.RowVector(i,0,n_parameters-1);
      m+=v;
      variance+=OuterProduct(v,v);
    }
  m=(1.0/in.rows)*m;
  variance=(1.0/in.rows)*variance - OuterProduct(m,m);

  // Compute decomposition
  Eig(d,U,0.5*(variance + Transpose(variance)));
  for (int j=d.dim-1; j >= 0; j--)  d.vector[j]=sqrt(d.vector[j]);
}

void TStorage::ComputeMean(TDenseVector &m)
{
  // Compute mean
  TDenseVector v;
  m=in.RowVector(0,0,n_parameters-1);
  for (int i=in.rows-1; i >= 0; i--)
    {
      v=in.RowVector(i,0,n_parameters-1);
      m+=v;
    }
  m=(1.0/in.rows)*m;
}

void TStorage::SetupForInput(int level, int number_rings)
{
  // set info
  level_in=level;
  ring_offset.resize(number_rings);
  ring_width.resize(number_rings);
  ring_boundary.Resize(number_rings-1);

  // read file
  fstream in_file;
  string filename=OpenFile(in_file,false,"draws",level);
  try
    {
      int count;
      in_file >> count;
      in.UniqueMemory(count,size,false);
      in_file >> in;
    }
  catch(...)
    {
      cerr << "SetupForInput(): Error reading " << filename << endl;
      abort();
    }
  in_file.close();
  in.SortRows(n_parameters);

  for (int m=1+(in.rows)/number_rings, hi=in.rows, i=number_rings-1; i >= 0; i--)
    {
      if (hi == 0)
	{
	  cerr << "SetupForInput(): No draws in some rings" << endl;
	  abort();
	}
      int j = (hi > m) ? hi-m : 0;
      for ( ; (j > 0) && (in(j,n_parameters) == in(j-1,n_parameters)); j--);
      ring_width[i]=hi-j;
      ring_offset[i]=hi=j;
    }

  for (int i=0; i < ring_boundary.dim; i++)
    ring_boundary[i]=in(ring_offset[i+1],n_parameters);
}

void TStorage::SetupForInput(int level, int number_rings, double temperature)
{
  SetupForInput(level,number_rings);
  K_level_in_minus_one=temperature;
  SetupWeights();
}

double TStorage::SetupForInput(int level, int number_rings, double K_max, double min_ess)
{
  SetupForInput(level,number_rings);
  SetTemperature(K_max, min_ess);
  SetupWeights();
  return K_level_in_minus_one;
}

void TStorage::EqualWeightedDraw(TDraw &x, int ring)
{
  if ((ring < 0) || (ring >= (int)ring_offset.size()))
    {
      cerr << "EqualWeightedDraw(): invalid ring" << endl;
      abort();
    }
  int idx=ring_offset[ring] + (int)floor(ring_width[ring] * dw_uniform_rnd());
  if (idx >= ring_offset[ring]+ring_width[ring]) idx=ring_offset[ring]+ring_width[ring]-1;
try
  {
    x.Set(in.RowVector(idx,0,n_parameters-1), in(idx,n_parameters), in(idx,n_parameters+1), in(idx,n_parameters+2), in(idx,n_parameters+3));

    // test ring
    if ((ring != Ring(x.LogPosterior())) || (x.LogPosterior() < LowerRingBoundary(ring)) || (UpperRingBoundary(ring) <= x.LogPosterior()))
      {
	cerr << setprecision(9) << scientific;
	cerr << LowerRingBoundary(ring) << " <= " << x.LogPosterior() << " <= " << UpperRingBoundary(ring) << " violated?  desired ring = "
	     << ring << "  computed ring = " << Ring(x.LogPosterior()) << endl;
	cerr << setprecision(9) << scientific << x.LogPosterior() - LowerRingBoundary(ring) << "   " << UpperRingBoundary(ring) - x.LogPosterior()
	     << " both should be positive" << endl;
	cerr << in(idx-1,n_parameters) << "  " << in(idx,n_parameters) << "  " << in(idx+1,n_parameters) << endl;
	abort(); 
      }
  }
 catch(dw_exception &e)
   {
     cerr << "EqualWeightedDraws(): exception - " << e.what() << endl;
     cerr << "idx " << idx << " - n_parameters " << n_parameters << " - rows " << in.rows << " - cols " << in.cols << endl;
     throw;
   }
}

void TStorage::ImportanceWeightedDraw(TDraw &x)
{
  int idx = UpperBound(cumulative_importance_weights, dw_uniform_rnd());
  if (idx >= cumulative_importance_weights.dim) idx=cumulative_importance_weights.dim;
  x.Set(in.RowVector(idx,0,n_parameters-1), in(idx,n_parameters), in(idx,n_parameters+1), in(idx,n_parameters+2), in(idx,n_parameters+3));
}

void TStorage::FreeInput(void)
{
  in.Resize(0,size);
  importance_weights.Resize(0);
  cumulative_importance_weights.Resize(0);
  ring_offset.resize(0);
  ring_width.resize(0);
  ring_boundary.Resize(0);
  level_in=-1;
}

void TStorage::SetupForOutput(int level, int node)
{
  level_out=level;
  node_out=node;
  idx_out=0;
  file_number_out=0;
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

void TStorage::WriteDraws(void)
{
  fstream out_file;
  string filename=OpenFile(out_file,true,"draws",level_out,node_out,file_number_out);
  try
    {
      out_file << setprecision(9) << scientific;
      out_file << idx_out << endl;
      out_file << out;
      idx_out=0;
      file_number_out++;
    }
  catch(...)
    {
      cerr << "Unable to write output to " << filename << endl;
      abort();
    }
  out_file.close();
}

void TStorage::ConsolidateDraws(int level_in, int level_out)
{
  string filename;
  fstream file;
  vector<TDenseMatrix> draws;
  int total_count=0;
  int node=0;
  while (1)
    {
      int file_number=-1;
      while (1)
	{
	  filename=MakeFullFileName("draws",level_in,node,file_number+1);
	  file.open(filename.c_str(), ios::in);
	  if (!file.is_open()) break;
	  try
	    {
	      int count;
	      file >> count;
	      TDenseMatrix draw(count,size);
	      file >> draw;
	      draws.push_back(draw);
	      total_count+=count;
	      file.close();
	      remove(filename.c_str()); 
	    }
	  catch(...)
	    {
	      cerr << "Error reading " << filename << endl;
	      abort();
	    }
	  file_number++;
	}
      if (file_number == -1) break;
      node++;
    }

  if (node == 0)
    {
      cerr << "No draws files to consolidate" << endl;
      abort();
    }

  OpenFile(file,true,"draws",level_out);
  try
    {
      file << setprecision(9) << scientific;
      file << total_count << endl;
      for (int ii=0; ii < draws.size(); ii++)
	file << draws[ii];
    }
  catch(...)
    {
      cerr << "Error writing to " << filename << endl;
      abort();
    }
  file.close();
}
/////////////////////////////////////////////////////////////////////////////////
/// End TStorage
/////////////////////////////////////////////////////////////////////////////////

void TStorage::PrintWeights(void)
{
  // check weights
  for (int i=0; i < importance_weights.dim; i++)
    {
      if ((importance_weights[i] < 0) || (cumulative_importance_weights[i] < 0))
	{
	  cerr << "negative weights" << endl;
	  abort();
	}
 
      if (i > 0) 
	{
	  if (fabs(cumulative_importance_weights[i] - cumulative_importance_weights[i-1] - importance_weights[i]) > 1.0e-7)
	    {
	      cerr << "cumulative weights not the sum of weights" << endl;
	      abort();
	    }
	}
      else if (fabs(cumulative_importance_weights[i] - importance_weights[i]) > 1.0e-7)
	{
	  cerr << "cumulative weights not the sum of weights" << endl;
	  abort();
	}
    }

  fstream file;
  OpenFile(file,true,"weights",level_in);
  for (int i=0; i < importance_weights.dim; i++) file << importance_weights[i] << " " << cumulative_importance_weights[i] << endl;
  file.close();

  OpenFile(file,true,"ring_info",level_in);
  for (int i=0; i < ring_offset.size(); i++) file << ring_offset[i] << " "; file << endl << endl;
  for (int i=0; i < ring_width.size(); i++) file << ring_width[i] << " "; file << endl << endl;
  for (int i=0; i < ring_boundary.dim; i++) file << in(ring_offset[i]+ring_width[i],n_parameters) << " "; file << endl << endl;
  for (int i=0; i < ring_boundary.dim; i++) file << ring_boundary[i] << " "; file << endl << endl;
  file.close();
}
