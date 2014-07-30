#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "CSampleIDWeight.h"
#include "CPutGetBin.h"
#include "CStorageHead.h"

using namespace std;

CStorageHead::CStorageHead(int size_each_data, int _node_index, const string & _run_id, size_t _storage_marker, string _file_location, size_t _number_level) : cluster_node(_node_index), run_id(_run_id), storage_marker(_storage_marker), filename_base(_file_location), bin(vector<vector<CPutGetBin> >(_number_level+1)), energy_lower_bound(vector<vector<double> >(_number_level+1)) 
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
		bin[i].push_back(CPutGetBin(size_each_data, bin_id_string.str(),0,storage_marker,filename_base+str.str(), cluster_node)); 
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

void CStorageHead::ClearSample(int level)
{
	if (level >= (int)bin.size())
	{
		cerr << "CStorageHead::ClearSample() : level index exceeds the range.\n"; 
		abort(); 
	}
	for (int j=0; j<(int)bin[level].size(); j++)
		bin[level][j].DisregardHistorySamples(); 
}

bool compare_CSampleIDWeight_BasedOnEnergy(const CSampleIDWeight &i, const CSampleIDWeight &j)
{
	return -i.weight < -j.weight; 
}

size_t CStorageHead::binning_equal_size(int level, size_t bin_number)
{
	if (level >=(int)bin.size())
        {
                cerr << "CStorageHead::binning_equal_size() : level index exceeds the range.\n";
                abort();
        }
	if (bin[level].size() > 1)
        {
                cerr << "CStorageHead::binning_equal_size() : it seems that binning has been done before.\n";
                abort();
        }	
        vector<CSampleIDWeight> sample;
        if (!DrawAllSample(level, sample) )
        {
                cerr << "CStorageHead::binning_equal_size() : error occurred when loading all samples.\n";
                abort();
        }
	sort(sample.begin(), sample.end(), compare_CSampleIDWeight_BasedOnEnergy); 
	// ascending energy (equivalently, descending on weight)
	int size_each_data = (int)sample[0].GetSize_Data(); 
        DisregardHistorySamples(level);

        stringstream str, bin_id_string;
        str << run_id << "/" << run_id << ".binary/";
        bin[level].clear();
	
	size_t sBin = (size_t)ceil((double)(sample.size())/(double)bin_number);
	
	int iSample=0, iBin=0;
	while (iSample < (int)(sample.size())) 
	{
		energy_lower_bound[level].push_back(-sample[iSample].weight); 	
		bin_id_string.str(string());
                bin_id_string << level << "." << iBin;
                bin[level].push_back(CPutGetBin(size_each_data, bin_id_string.str(),0,storage_marker,filename_base+str.str(), cluster_node));

		sBin = iSample+sBin == (int)(sample.size())-1 ? sBin+1 : sBin;
                for (int j=iSample; j<iSample+(int)sBin && j<(int)(sample.size()); j++)
                	bin[level][iBin].DepositSample(sample[j]);
               	iSample += sBin;
               	iBin++;
	}
	sample.clear(); 
	return bin[level].size(); 
}

bool CStorageHead::DrawLeastWeightSample(int level, int _bin_id, CSampleIDWeight &sample) 
{
	if (level >= (int)bin.size() || _bin_id >=(int)bin[level].size() )
	{
		cerr << "CStorageHead::DrawLeastWeightSample : level or bin index exceeds the range.\n"; 
		abort(); 
	}
	return bin[level][_bin_id].DrawLeastWeightSample(sample); 
}

bool CStorageHead::Draw_K_LeastWeightSample(size_t K, int level, int _bin_id, vector<CSampleIDWeight> &sample)
{
	if (level >= (int)bin.size() || _bin_id >= (int)bin[level].size() )
	{
		cerr << "CStorageHead::Draw_K_LeastWeightSample : level or bin index exceeds the range.\n"; 
		abort(); 
	}
        return bin[level][_bin_id].Draw_K_LeastWeightSample(K, sample);
}

bool CStorageHead::DrawMostWeightSample(int level, int _bin_id, CSampleIDWeight &sample) 
{
	if (level >=(int)bin.size() || _bin_id >=(int)bin[level].size() )
	{
		cerr << "CStorageHead::DrawMostWeightSample : level or bin index exceeds the range.\n"; 
		abort(); 
	}
	return bin[level][_bin_id].DrawMostWeightSample(sample); 
}

bool CStorageHead::Draw_K_MostWeightSample(size_t K, int level, int _bin_id, vector<CSampleIDWeight> &sample)
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

bool CStorageHead::DrawAllSample(int level, vector<CSampleIDWeight>&sample)
{
	if (level >= (int)bin.size())
	{
		cerr << "CStorageHead::DrawAllSample : level index exceeds the range.\n"; 
		abort(); 
	}
	for (int ii=0; ii<(int)(bin[level].size()); ii++)
	{
		if (!bin[level][ii].DrawAllSample(sample))		
			return false; 
	}
	return true; 
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

void CStorageHead::ClearStatus(int level)
{
	for (int iBin=0; iBin<(int)(bin[level].size()); iBin++)
		bin[level][iBin].ClearStatus(); 
}
