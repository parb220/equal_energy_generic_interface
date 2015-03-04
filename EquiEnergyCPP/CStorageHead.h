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
	int cluster_node;
	string run_id; 
	size_t storage_marker; 
	string filename_base; 
	
	vector<vector<CPutGetBin> > bin;
	vector<vector<double> > energy_lower_bound; 
public: 
	CStorageHead(int _node_index=0, const string & _run_id=string(), size_t _storage_marker=10000, string _file_location=string(), size_t _number_stage=1); 
	~CStorageHead(); 

	bool makedir(); 
	
	int BinIndex(int stage, double energy) const;  
	int DepositSample(int stage, int _bin_id, const CSampleIDWeight &) ;
	void finalize(int stage);
	void consolidate(int stage); 
	void restore(int stage);
	void ClearDepositDrawHistory(int stage); 
	void ClearSample(int stage); 

	size_t binning_equal_size(int stage, size_t bin_number, bool if_unstructured=false, int data_size=1); 

	bool DrawLeastWeightSample(int stage, int bin_id, CSampleIDWeight &); 
	bool Draw_K_LeastWeightSample(size_t, int stage, int bin_id, vector<CSampleIDWeight> &) ; 
	bool DrawMostWeightSample(int stage, int bin_id, CSampleIDWeight &); 
	bool Draw_K_MostWeightSample(size_t, int stage, int bin_id, vector<CSampleIDWeight> &) ;
	bool DrawSample(int stage, int bin_id, CSampleIDWeight &) ; 
	bool DrawAllSample(int stage, vector<CSampleIDWeight> &, bool unstructured=false, int data_size=1) ; 
	
	size_t GetNumberRecrod(int stage, int index) const ;  

	/* for reassigning samples into different bins */
	void DisregardHistorySamples(int stage); 
	void RestoreForFetch(int stage); 
	bool empty(int stage) const; 
	
	/* reset setting */
	void ClearStatus(int ); 

	/* get set energy lower bound */
 	int GetNumber_Bin(int stage) const;
        double GetEnergyLowerBound(int stage, int index) const;
        void ResizeBin(int stage, int number, int data_size); 
        void SetEnergyLowerBound(int stage, int index, double e); 
	void ClearBin(int stage); 
	void InitializeBin(int stage, int data_size); 

}; 

#endif
