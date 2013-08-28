#ifndef _STORAGE_HEAD_
#define _STORAGE_HEAD_

#include <vector>
#include <string>

using namespace std; 

class CSampleIDWeight; 
class CPutGetBin; 

class CStorageHead
{
protected:
	unsigned int cluster_node;
	unsigned int run_id; 
	size_t storage_marker; 
	string filename_base; 
	
	vector<vector<CPutGetBin> > bin;
	vector<vector<double> > energy_lower_bound; 
public: 
	CStorageHead(unsigned int _node_index=0, unsigned int _run_id=0, size_t _storage_marker=10000, string _file_location=string(), size_t _number_level=1); 
	~CStorageHead(); 

	bool makedir(); 
	size_t Number_Bin(unsigned int level) const; 
	
	unsigned int BinIndex(unsigned int level, double energy) const;  
	unsigned int DepositSample(unsigned int level, unsigned int _bin_id, const CSampleIDWeight &) ;
	void finalize(unsigned int level);
	void consolidate(unsigned int level); 
	void restore(unsigned int level);
	void ClearDepositDrawHistory(unsigned int level); 

	size_t binning(unsigned int level, size_t bin_number_lb, double bin_width_ub); 

	bool DrawLeastWeightSample(unsigned int level, unsigned int bin_id, CSampleIDWeight &) const; 
	bool Draw_K_LeastWeightSample(size_t, unsigned int level, unsigned bin_id, vector<CSampleIDWeight> &) const; 
	bool DrawMostWeightSample(unsigned int level, unsigned int bin_id, CSampleIDWeight &)const; 
	bool Draw_K_MostWeightSample(size_t, unsigned int level, unsigned int bin_id, vector<CSampleIDWeight> &) const;
	bool DrawSample(unsigned int level, unsigned bin_id, CSampleIDWeight &) ; 
	
	size_t GetNumberRecrod(unsigned int level, unsigned int index) const ;  

	/* for reassigning samples into different bins */
	virtual void DisregardHistorySamples(unsigned int level); 
	void RestoreForFetch(unsigned int level); 
	bool empty(unsigned int level) const; 
}; 

#endif
