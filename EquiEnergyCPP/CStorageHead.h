#ifndef _STORAGE_HEAD_
#define _STORAGE_HEAD_

#include <vector>
#include <string>

using namespace std; 

class CSampleIDWeight; 
class CPutGetBin; 
bool Solve_Polynomial_Equation(vector<double> &, size_t n, double t0, double tk_1); 

class CStorageHead
{
protected:
	int cluster_node;
	string run_id; 
	size_t storage_marker; 
	string filename_base; 
	
	vector<vector<CPutGetBin> > bin;
	vector<vector<double> > energy_lower_bound; 
public: 
	CStorageHead(int size_each_data, int _node_index=0, const string & _run_id=string(), size_t _storage_marker=10000, string _file_location=string(), size_t _number_level=1); 
	~CStorageHead(); 

	bool makedir(); 
	size_t Number_Bin(int level) const; 
	
	int BinIndex(int level, double energy) const;  
	int DepositSample(int level, int _bin_id, const CSampleIDWeight &) ;
	void finalize(int level);
	void consolidate(int level); 
	void restore(int level);
	void ClearDepositDrawHistory(int level); 
	void ClearSample(int level); 

	size_t binning_equal_size(int level, size_t bin_number); 

	bool DrawLeastWeightSample(int level, int bin_id, CSampleIDWeight &); 
	bool Draw_K_LeastWeightSample(size_t, int level, int bin_id, vector<CSampleIDWeight> &) ; 
	bool DrawMostWeightSample(int level, int bin_id, CSampleIDWeight &); 
	bool Draw_K_MostWeightSample(size_t, int level, int bin_id, vector<CSampleIDWeight> &) ;
	bool DrawSample(int level, int bin_id, CSampleIDWeight &) ; 
	bool DrawAllSample(int level, vector<CSampleIDWeight> &, bool unstructured=false, int data_size=1) ; 
	
	size_t GetNumberRecrod(int level, int index) const ;  

	/* for reassigning samples into different bins */
	void DisregardHistorySamples(int level); 
	void RestoreForFetch(int level); 
	bool empty(int level) const; 
	
	/* reset setting */
	void ClearStatus(int ); 
}; 

#endif
