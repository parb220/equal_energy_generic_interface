#ifndef _CSAMPLE_ID_WEIGHT
#define _CSAMPLE_ID_WEIGHT

#include <iostream>
#include <fstream>
#include <vector>
#include "dw_dense_matrix.hpp"

using namespace std; 

class CEquiEnergyModel; 

class CSampleIDWeight
{
protected:
	bool calculated; 
public:
	TDenseVector data; 
	int id; 
	double weight; 

	CSampleIDWeight(const TDenseVector &_x=TDenseVector(), int _id=0, double _weight=0.0, bool _calculated=false); 
	CSampleIDWeight(const CSampleIDWeight &); 
	~CSampleIDWeight();

	CSampleIDWeight & operator=(const CSampleIDWeight &);
	bool PartialCopyFrom(const CSampleIDWeight &, unsigned int offset, size_t length);
	bool PartialCopyFrom(unsigned int offset1, const CSampleIDWeight &, unsigned int offset2, size_t length);

	size_t GetSize_Data() const; 
	void DataChanged(); 

	friend istream& read(istream &, CSampleIDWeight *); 
	friend ostream& write(ostream &, const CSampleIDWeight *);

	friend istream& operator>>(istream &, CSampleIDWeight &); 
	friend ostream& operator<<(ostream &, const CSampleIDWeight &); 
	friend class CEquiEnergyModel; 
}; 

bool LoadSampleFromFile(const string &file, vector<CSampleIDWeight> &Y); 
bool SaveSampleToFile(const string &file, const vector<CSampleIDWeight> &Y); 

#endif
