#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <cstdlib>
#include <cmath>

#include "CSampleIDWeight.h"
#include "CPutGetBin.h"
#include "CStorageHead.h"

using namespace std;

CStorageHead::CStorageHead(int _node_index, const string & _run_id, size_t _storage_marker, string _file_location, size_t _number_level) : cluster_node(_node_index), run_id(_run_id), storage_marker(_storage_marker), filename_base(_file_location), bin(vector<vector<CPutGetBin> >(_number_level+1)), energy_lower_bound(vector<vector<double> >(_number_level+1)) 
{
	stringstream str; 
	str << run_id << "/" << run_id << ".binary/";

	// bin[i]: bins for the i-th temperature level
	// bin[i][j]: j-th bin for the i-th temperature level
	// energy_lower_bound[i]: energy lower bounds for the i-th temperature level
	// energy_lower_bound[i][j]: energy lower bound for the j-th bin of the i-th level
	
	stringstream bin_id_string; 
	// Each level has at least one bin for depositing
	for (int i=0; i<(int)bin.size(); i++)
	{
		bin_id_string.str(string()); 
		bin_id_string << i << ".0"; 
		bin[i].push_back(CPutGetBin(bin_id_string.str(),0,storage_marker,filename_base+str.str(), cluster_node)); 
	}
}

CStorageHead::~CStorageHead()
{
}

bool CStorageHead::makedir()
{
	stringstream str; 
	str.str(string()); 
	str << run_id; 
	string dir = filename_base + str.str(); 
	errno = 0; // clear error number
	int status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (status !=0 && errno != EEXIST)
		return false;
 
	str.str(string()); 
	str << run_id << "/" << run_id << ".binary" ; 
	dir = filename_base + str.str();
	errno = 0; 
	status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); 
	if (status == 0 || errno == EEXIST)
		return true; 
	else 
		return false; 
}

size_t CStorageHead::Number_Bin(int level) const
{
	if (level >=(int) bin.size())
	{
		cerr << "CStorageHead::Number_Bin : level index exceeds the range.\n"; 
		abort(); 
	}
	return bin[level].size(); 
}

int CStorageHead::BinIndex(int level, double energy) const 
{
	if (level >= (int)energy_lower_bound.size())
	{
		cerr << "CStorageHead::energy_lower_bound.resize(bin_number); BinIndex() : level index exceeds the range.\n"; 
		abort(); 
	}
	if (energy_lower_bound[level].empty())
		return 0; 
	for (int j=1; j<(int)energy_lower_bound[level].size(); j++)
	{
		if (energy < energy_lower_bound[level][j])
			return j-1; 
	}
	return energy_lower_bound[level].size()-1; 
}

int CStorageHead::DepositSample(int level, int _bin_id, const CSampleIDWeight &sample)
{
	if (level >= (int)bin.size())
	{
		cerr << "CStorageHead::DepositSample() : level index exceeds the range.\n"; 
		abort(); 
	}
	return bin[level][_bin_id].DepositSample(sample); 
}

void CStorageHead::finalize(int level)
{
	if (level >= (int)bin.size())
	{
		cerr << "CStorageHead::finalize() : level index exceeds the range.\n"; 
		abort(); 
	}
	for (int j=0; j<(int)bin[level].size(); j++)
		bin[level][j].finalize(); 
}

void CStorageHead::consolidate(int level)
{
	if (level >=(int)bin.size())
	{
		cerr << "CStorageHead::consolidate() : level index exceeds the range.\n"; 
		abort(); 
	}
	for (int j=0; j<(int)bin[level].size(); j++)
		bin[level][j].consolidate(); 
}

void CStorageHead::restore(int level)
{
	if (level >= (int)bin.size())
	{
		cerr << "CStorageHead::restore() : level index exceeds the range.\n"; 
                abort();
	}
	for (int j=0; j<(int)bin[level].size(); j++)
		bin[level][j].restore(); 
}

void CStorageHead::ClearDepositDrawHistory(int level)
{
	if (level >= (int)bin.size())
	{
		cerr << "CStorageHead::ClearDepositDrawHistory() : level index exceeds the range.\n"; 
		abort(); 
	}
	for (int j=0; j<(int)bin[level].size(); j++) 
		bin[level][j].ClearDepositDrawHistory(); 
}

size_t CStorageHead::binning_geometric(int level, size_t bin_number)
{
	if (level >=(int)bin.size())
	{
		cerr << "CStorageHead::binning_geometric() : level index exceeds the range.\n"; 
		abort(); 
	}
	if (bin[level].size() > 1)
	{
		cerr << "CStorageHead::binning_geometric() : it seems that binning has been done before.\n";
		abort(); 
	}
	size_t nSample=bin[level][0].GetTotalNumberRecord();
        vector<CSampleIDWeight> sample;
        if (!Draw_K_MostWeightSample(nSample, level, 0, sample) )
        {
                cerr << "CStorageHead::binning_geometric() : error occurred when loading all samples.\n";
                abort();
        }
        DisregardHistorySamples(level);
	
	stringstream str, bin_id_string;
        str << run_id << "/" << run_id << ".binary/";
	bin[level].clear();
	for (int i=0; i<(int)bin_number; i++)
	{
		bin_id_string.str(string());
                bin_id_string << level << "." << i;
                bin[level].push_back(CPutGetBin(bin_id_string.str(),0,storage_marker,filename_base+str.str(), cluster_node));
	}
	double energy_min = -sample[0].weight, energy_max=-sample[nSample-1].weight;
	// double energy_min = 1950, energy_max =1975; 
	if (!Solve_Polynomial_Equation(energy_lower_bound[level], bin_number,energy_min,energy_max) ) 
	{
		cerr << "CStorageHead::binning_geometric() : polynomial equation cannot be solved.\n"; 
		abort(); 
	}
	for(int i=0; i<(int)nSample; i++)
		DepositSample(level, BinIndex(level,-sample[i].weight), sample[i]); 
	return bin_number; 
}

size_t CStorageHead::binning(int level, size_t bin_number_lb, double bin_width_ub)
{
	// (1) Before binning is carried out, bin[level] has only one bin containing all
	// samples that have been drawn so far. 
	if (level >= (int) bin.size())
	{
		cerr << "CStorageHead::binning() : level index exceeds the range.\n"; 
		abort(); 
	}
	if (bin_width_ub <= 0)
	{
		cerr << "CStorageHead::binning() : upper-bound for the bin width has to be positive.\n"; 
		abort(); 
	}
	if (bin[level].size() > 1)
	{
		cerr << "CStorageHead::binning() : it seems that binning has been done before.\n"; 
		abort(); 
	}
	size_t nSample=bin[level][0].GetTotalNumberRecord(); 
	vector<CSampleIDWeight> sample; 
	if (!Draw_K_MostWeightSample(nSample, level, 0, sample) )
	{
		cerr << "CStorageHead::binning() : error occurred when loading all samples.\n"; 
		abort(); 
	}
	DisregardHistorySamples(level); 

	stringstream str; 
	str << run_id << "/" << run_id << ".binary/";
	// number of bins has to be at least number of levels
	size_t nBin = bin_number_lb > bin.size() ? bin_number_lb : bin.size(); 
	// initial guess of the size of bins is based on nBin
	size_t sBin_initial = (size_t)ceil((double)(nSample)/(double)(nBin));
	
	bin[level].clear(); // Get rid of bin[0] because all samples are now in sample
	energy_lower_bound[level].clear(); 
	int iSample=0, iBin=0; 
	stringstream bin_id_string; 
	while (iSample < (int)nSample)
	{
		size_t sBin = iSample+sBin_initial<nSample ? sBin_initial : nSample-iSample-1; 
		// width of a bin = sample[iSample].weight - sample[iSample+sBin].weight
		// bin width has to be less than or equal to bin_width_ub
		/*while (sBin>1 && fabs(sample[iSample].weight-sample[iSample+sBin].weight) > bin_width_ub)
			sBin = sBin/2; 
		*/
	
		// Now the bin boundary has been determined. 
		// It is possible that sBin=1 ???????
		energy_lower_bound[level].push_back(-sample[iSample].weight); 
		bin_id_string.str(string()); 
		bin_id_string << level << "." << iBin; 
		bin[level].push_back(CPutGetBin(bin_id_string.str(),0,storage_marker,filename_base+str.str(), cluster_node)); 
	
		// If this is the last bin, then it has to include the last element of sample	
		sBin = (iSample+sBin == nSample-1) ? sBin+1 : sBin; 
		for (int j=iSample; j<iSample+(int)sBin; j++)
			bin[level][iBin].DepositSample(sample[iSample]); 

		iSample += sBin; 
		iBin++; 	
	}
	return bin[level].size(); 
}

bool CStorageHead::DrawLeastWeightSample(int level, int _bin_id, CSampleIDWeight &sample) const
{
	if (level >= (int)bin.size() || _bin_id >=(int)bin[level].size() )
	{
		cerr << "CStorageHead::DrawLeastWeightSample : level or bin index exceeds the range.\n"; 
		abort(); 
	}
	return bin[level][_bin_id].DrawLeastWeightSample(sample); 
}

bool CStorageHead::Draw_K_LeastWeightSample(size_t K, int level, int _bin_id, vector<CSampleIDWeight> &sample) const
{
	if (level >= (int)bin.size() || _bin_id >= (int)bin[level].size() )
	{
		cerr << "CStorageHead::Draw_K_LeastWeightSample : level or bin index exceeds the range.\n"; 
		abort(); 
	}
        return bin[level][_bin_id].Draw_K_LeastWeightSample(K, sample);
}

bool CStorageHead::DrawMostWeightSample(int level, int _bin_id, CSampleIDWeight &sample) const
{
	if (level >=(int)bin.size() || _bin_id >=(int)bin[level].size() )
	{
		cerr << "CStorageHead::DrawMostWeightSample : level or bin index exceeds the range.\n"; 
		abort(); 
	}
	return bin[level][_bin_id].DrawMostWeightSample(sample); 
}

bool CStorageHead::Draw_K_MostWeightSample(size_t K, int level, int _bin_id, vector<CSampleIDWeight> &sample) const
{
	if (level >=(int)bin.size() || _bin_id >= (int)bin[level].size() )
        {
                cerr << "CStorageHead::Draw_K_MostWeightSample : level or bin index exceeds the range.\n";
                abort();
        }

        return bin[level][_bin_id].Draw_K_MostWeightSample(K, sample);
}

bool CStorageHead::DrawSample(int level, int _bin_id, CSampleIDWeight &sample)
{
	if (level >= (int)bin.size() || _bin_id >= (int)bin[level].size() )
        {
                cerr << "CStorageHead::DrawSample : level or bin index exceeds the range.\n";
                abort();
        }

	return bin[level][_bin_id].DrawSample(sample); 	
}

size_t CStorageHead::GetNumberRecrod(int level, int index) const
{
	if (level >= (int)bin.size() || index >= (int)bin[level].size() )
	{
		cerr << "CStorageHead::GetNumberRecrod : level or bin index exceeds the range.\n";
                abort();
	}
	return bin[level][index].GetTotalNumberRecord();
}

void CStorageHead::DisregardHistorySamples(int level)
{
	if (level >= (int)bin.size())
	{
		cerr << "CStorageHead::DisregardHistorySamples : level index exceeds the range.\n"; 
		abort(); 
	}
	for (int j=0; j<(int)bin[level].size(); j++)
		bin[level][j].DisregardHistorySamples();
}

bool CStorageHead::empty(int level) const
{
	if (level >= (int)bin.size() )
	{
		cerr << "CStorageHead::empty : level index exceeds the range.\n"; 
		abort(); 
	}
	for (int j=0; j<(int)bin[level].size(); j++)
	{
		if (!bin[level][j].empty())
			return false; 
	}
	return true; 
}

void CStorageHead::RestoreForFetch(int level)
{
	if (level >=(int) bin.size() )
	{
		cerr << "CStorageHead::RestoreForFetch : level index exceeds the range.\n";
		abort(); 
	}
	for (int j=0; j<(int)bin[level].size(); j++)
		bin[level][j].RestoreForFetch(); 
}

