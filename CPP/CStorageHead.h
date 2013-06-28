#ifndef _STORAGE_HEAD_
#define _STORAGE_HEAD_

#include <vector>
#include <string>
#include "CPutGetBin.h"
#include "CSampleIDWeight.h"

using namespace std; 

class CStorageHead
{
protected:
	int cluster_node;
	int run_id; 
	size_t storage_marker; 
	size_t number_bins; 
	vector <CPutGetBin> bin; 
	string filename_base; 
public: 
	CStorageHead(int _run_id = 0, size_t _storage_marker = 10000, size_t _number_bins = 1, string file=string(), int _node_index=0); 
	virtual ~CStorageHead(); 

	virtual unsigned int DepositSample(unsigned int, const CSampleIDWeight &) ; 
	virtual bool DrawLeastWeightSample(unsigned int, CSampleIDWeight &) const; 
	virtual bool DrawMostWeightSample(unsigned int, CSampleIDWeight &)const; 
	virtual bool DrawSample(unsigned int, CSampleIDWeight &) ; 
	bool makedir(); 
	void finalize(unsigned int =1,  unsigned int =0)  ;
	size_t GetNumberBins() const { return number_bins; }
	size_t GetNumberRecrod(unsigned int index) const { return bin[index].GetTotalNumberRecord(); }

	/* for reassigning samples into different bins */
	virtual void DisregardHistorySamples(unsigned int = 1, unsigned int = 0); 
	virtual void ClearDepositDrawHistory(unsigned int = 1, unsigned int = 0); 
	virtual void consolidate(unsigned int = 1, unsigned int = 0); 
	void restore(unsigned int = 1, unsigned int = 0);
	void RestoreForFetch(unsigned int = 1, unsigned int = 0); 
	bool empty(unsigned int = 1, unsigned int = 0) const; 
}; 

#endif
