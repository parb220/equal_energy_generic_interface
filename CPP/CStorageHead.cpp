#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "CStorageHead.h"

using namespace std;

CStorageHead::CStorageHead(int _run_id, size_t _storage_marker, size_t _number_bins, string file_location, int _node_index) : cluster_node(_node_index), run_id(_run_id), storage_marker(_storage_marker), number_bins(_number_bins), filename_base(file_location)
{
	stringstream str; 
	str << run_id << "/" << run_id << ".binary/";
	for (unsigned int i=0; i<bin.size(); i++)
		bin[i] = CPutGetBin(0,0,storage_marker,filename_base+str.str(), cluster_node); 
}

CStorageHead::~CStorageHead()
{
}


unsigned int CStorageHead::DepositSample(unsigned int _bin_id, const CSampleIDWeight &sample)
{
	return bin[_bin_id].DepositSample(sample); 
}


bool CStorageHead::DrawLeastWeightSample(unsigned int _bin_id, CSampleIDWeight &sample) const
{
	return bin[_bin_id].DrawLeastWeightSample(sample); 
}

bool CStorageHead::DrawMostWeightSample(unsigned int _bin_id, CSampleIDWeight &sample) const
{
	return bin[_bin_id].DrawMostWeightSample(sample); 
}

bool CStorageHead::DrawSample(unsigned int _bin_id, CSampleIDWeight &sample)
{
	return bin[_bin_id].DrawSample(sample); 	
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

void CStorageHead::finalize(unsigned int start_bin, unsigned int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (unsigned int i=start_bin; i<=end_bin; i++)
			bin[i].finalize(); 
	}
	else
	{
		for (unsigned int i=0; i<number_bins; i++)
			bin[i].finalize(); 
	}
}

void CStorageHead::DisregardHistorySamples(unsigned int start_bin, unsigned int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (unsigned int i=start_bin; i<=end_bin; i++)
			bin[i].DisregardHistorySamples();
	}
	else 
	{
		for (unsigned int i=0; i<number_bins; i++)
			bin[i].DisregardHistorySamples(); 
	} 
}

void CStorageHead::ClearDepositDrawHistory(unsigned int start_bin, unsigned int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (unsigned int i=start_bin; i<=end_bin; i++)
			bin[i].ClearDepositDrawHistory(); 
	}
	else 
	{
		for (unsigned int i=0; i<number_bins; i++)
			bin[i].ClearDepositDrawHistory(); 
	}
}

void CStorageHead::consolidate(unsigned int start_bin, unsigned int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (unsigned int i=start_bin; i<=end_bin; i++)
			bin[i].consolidate(); 
	}
	else 
	{
		for (unsigned int i=0; i<number_bins; i++)
			bin[i].consolidate(); 
	}
}

void CStorageHead::restore(unsigned int start_bin, unsigned int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
        {
                for (unsigned int i=start_bin; i<=end_bin; i++)
                        bin[i].restore();
        }
	else 
	{
		for (unsigned int i=0; i<number_bins; i++)
			bin[i].restore(); 
	}
}

bool CStorageHead::empty(unsigned int start_bin, unsigned int end_bin) const
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (unsigned int i=end_bin; i>=start_bin; i--)
			if (!bin[i].empty())
				return false; 
		return true; 
	}
	else 
	{
		for (unsigned int i=number_bins-1; i>=0; i--)
			if (!bin[i].empty())
				return false; 
		return true; 
	}
}

void CStorageHead::RestoreForFetch(unsigned int start_bin, unsigned int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
        {
                for (unsigned int i=start_bin; i<=end_bin; i++)
                        bin[i].RestoreForFetch();
        }
	else 
	{
		for (unsigned int i=0; i<number_bins; i++)
			bin[i].RestoreForFetch(); 
	}

}
