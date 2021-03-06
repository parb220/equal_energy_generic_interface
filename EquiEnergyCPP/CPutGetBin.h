#ifndef _PUT_GET_BIN_
#define _PUT_GET_BIN_

#include <string>
#include <vector>

class CSampleIDWeight; 

using namespace std; 

class CPutGetBin					 
{
protected:
	int size_each_data;  
	int suffix;  
	string id; 
	int nDumpFile;			// total number of samples generated by far 
	size_t capacity; 			// capacity of this bin for put and get;
	int nPutUsed; 			// available space for put (will be reset as capacity after each dump, and will decrease as a sample is depoisited)
	int nGetUsed; 			// available data for get (will be reset as capacity after each fetch, and will decrease as a sample is drawn)	
	string filename_prefix; 		// run_id.id.indxx (index = nSampleGeneratedByFar/capacity); 
	vector <CSampleIDWeight> dataPut; 	// space for put
	vector <CSampleIDWeight> dataGet; 	// data for get 

	string GetFileNameForDump() const; 
	vector <string > GetFileNameForFetch() const; 
	vector <string > GetFileNameForConsolidate() const; 
	bool Dump(const string & =string()); 			// dump the current materials to file;
	bool Fetch(const vector<string> &); 
 	bool ReadFromOneFile(const string &, int &, const vector<int> &index);
	vector <CSampleIDWeight> ReadSampleFromFile(const string & ) const;  

	int NumberRecord(const string &) const; 
	bool LoadLeastWeightSample(const string &, CSampleIDWeight &) const; 
	bool LoadMostWeightSample(const string &, CSampleIDWeight &) const; 
	bool Load_K_LeastWeightSample(size_t, const string &, vector<CSampleIDWeight> &) const; 
	bool Load_K_MostWeightSample(size_t, const string &, vector<CSampleIDWeight> &) const; 
public: 
	CPutGetBin(int _size_each_data, const string & _id=string(), int _nDumpFile=0, size_t _capacity=0, string _grandPrefix=string(), int _suffix=0); 
	CPutGetBin(const CPutGetBin &); 
	CPutGetBin & operator=(const CPutGetBin &); 
	~CPutGetBin() ; 

	void SetBinID(string _id, int _suffix=0) { id = _id; suffix=_suffix; }
	const string & GetBinID() const { return id; }

	// void SetNumberSamplesGeneratedByFar(int _nTotalSamples) { nSamplesGeneratedByFar = _nTotalSamples; }
	// int GetNumberSamplesGeneratedByFar() const { return nDumpFile*capacityPut+nPutUsed; }

	void SetCapacity (size_t _capacity) { capacity = _capacity; }
	int GetCapacity() const { return capacity; }
	
	size_t GetNumberFileForFetch() const; 
	size_t GetNumberFileForDump() const; 
	size_t GetNumberFileForConsolidate() const; 
	
	void SetFileNamePrefix(const string &_grandPrefix) { filename_prefix = _grandPrefix; } 
	string GetFileNamePrefix() const { return filename_prefix;}	

	int DepositSample(const CSampleIDWeight &); 

	bool DrawLeastWeightSample(CSampleIDWeight &) const;  
	bool DrawMostWeightSample(CSampleIDWeight &) const; 
	bool Draw_K_LeastWeightSample(size_t, vector<CSampleIDWeight> &) const; 
	bool Draw_K_MostWeightSample(size_t, vector<CSampleIDWeight> &) const; 
	bool DrawSample(CSampleIDWeight &); 
	bool DrawAllSample(vector<CSampleIDWeight> &)const; 

	void finalize(); 	// save unsaved data
	void consolidate(); 	// conslidate partial sample files into complete sample files
	void restore();	// load data from a partial file
	void RestoreForFetch(); // load data from a partial file but will not update it later. This is used for single-thread mpi version so that for each level it will load partial files for its higher level for ee draw later
	bool empty() const; 

	/* for reassigning samples into different bins */ 
	void DisregardHistorySamples(); 
	void ClearDepositDrawHistory(); 

	/* to get the number of records in this bin */
	size_t GetTotalNumberRecord() const; 
}; 

#endif
