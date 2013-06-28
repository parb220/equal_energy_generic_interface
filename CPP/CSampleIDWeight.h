#ifndef _CSAMPLE_ID_WEIGHT
#define _CSAMPLE_ID_WEIGHT

#include <iostream>
#include <fstream>
#include "dw_dense_matrix.hpp"

using namespace std; 

class CSampleIDWeight
{
public:
	TDenseVector data; 
	int id; 
	double weight; 

	CSampleIDWeight(const TDenseVector &_x=TDenseVector(), int _id=0, double _weight=0.0); 
	CSampleIDWeight(const CSampleIDWeight &); 
	~CSampleIDWeight();

	CSampleIDWeight & operator=(const CSampleIDWeight &);
	bool PartialCopyFrom(const CSampleIDWeight &, int offset, int length);
	bool PartialCopyFrom(int offset1, const CSampleIDWeight &, int offset2, int length);

	unsigned int GetSize_Data()
	{
		return sizeof(int)+data.dim*sizeof(double)+sizeof(int)+sizeof(double);
		// data.dim, data.vector, id, weight 
	}

	friend istream& read(istream &, CSampleIDWeight *); 
	friend ostream& write(ostream &, const CSampleIDWeight *);

	friend istream& operator>>(istream &, CSampleIDWeight &); 
	friend ostream& operator<<(ostream &, const CSampleIDWeight &); 
}; 

#endif
