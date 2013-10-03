#include <string>
#include "dw_dense_matrix.hpp"
#include "CSampleIDWeight.h"

using namespace std; 

void CSampleIDWeight::DataChanged()
{
	calculated = false; 
}

size_t CSampleIDWeight::GetSize_Data() const
{
	return sizeof(int)+(size_t)data.dim*sizeof(double)+sizeof(int)+sizeof(double)+sizeof(bool);
	// data.dim, data.vector, id, weight, calculated
}

CSampleIDWeight::CSampleIDWeight(const TDenseVector &_x, int _id, double _weight, bool _calculated) : id(_id), weight(_weight), calculated(_calculated)
{
	data.CopyContent(_x); 
}

CSampleIDWeight::CSampleIDWeight(const CSampleIDWeight &right) : id(right.id), weight(right.weight), calculated(right.calculated)
{
	data.CopyContent(right.data); 
}

CSampleIDWeight::~CSampleIDWeight()
{
}

CSampleIDWeight &CSampleIDWeight::operator = (const CSampleIDWeight &right)
{
	data.CopyContent(right.data); 
	id = right.id; 
	weight = right.weight; 
	calculated = right.calculated; 
	return *this; 
}

bool CSampleIDWeight::PartialCopyFrom(const CSampleIDWeight &right, int offset, size_t length)
{
	if (data.dim < offset+(int)length || right.data.dim < offset+(int)length)
	{
		cerr << "CSampleIDWeight::PartialCopyFrom() : index exceeding dimension" << endl; 
		return false; 
	}

	for (int i=0; i<(int)length; i++)
		data[i+offset] = right.data[i+offset]; 	
	
	if (offset ==0 && (int)length == data.dim && (int)length == right.data.dim )
	{
		weight = right. weight; 
		calculated = right.calculated; 
	}

	return true; 
}

bool CSampleIDWeight::PartialCopyFrom(int offset1, const CSampleIDWeight &right, int offset2, size_t length)
{
	if (data.dim < offset1+(int)length || right.data.dim < offset2+(int)length )
	{
		cerr << "CSampleIDWeight::PartialCopyFrom() : index exceeding dimension" << endl; 
		return false; 
	}
	for (int i=0; i<(int)length; i++)
		data[offset1+i] = right.data[offset2+i]; 

	if (offset1 == 0 && offset2 == 0 && (int)length == data.dim && (int)length == right.data.dim)
	{
		weight = right.weight;
		calculated = right.calculated; 
	}
	
	return true; 
}

istream & read (istream & input_stream, CSampleIDWeight *x)
{
	int dim; 
	input_stream.read((char*)&(dim), sizeof(int)); 
	x->data.Resize(dim);
	double temp_x; 
	for (int i=0; i<dim; i++)
	{
		input_stream.read((char*)&temp_x, sizeof(double)); 
		x->data[i] = temp_x; 
	}
	input_stream.read((char*)&(x->id), sizeof(int)); 
	input_stream.read((char*)&(x->weight), sizeof(double)); 
	input_stream.read((char*)&(x->calculated), sizeof(bool)); 
	return input_stream; 
}

ostream& write(ostream & output_stream, const CSampleIDWeight *x)
{
	output_stream.write((char*)&(x->data.dim), sizeof(int)); 
	double temp_x; 
	for (int i=0; i<x->data.dim; i++)
	{
		temp_x = x->data[i]; 
		output_stream.write((char*)&temp_x, sizeof(double)); 
	}
	output_stream.write((char*)&(x->id), sizeof(int)); 
	output_stream.write((char*)&(x->weight), sizeof(double)); 
	output_stream.write((char*)&(x->calculated), sizeof(bool)); 
	return output_stream; 
}

istream& operator >> (istream &inputStr, CSampleIDWeight &sample)
{
	int dim; 
	inputStr >> dim; 
	sample.data.Resize(dim); 
	for (int i=0; i<sample.data.dim; i++)
		inputStr >> sample.data[i]; 
	inputStr >> sample.id; 
	inputStr >> sample.weight; 
	inputStr >> sample.calculated;
	return inputStr; 
}

ostream& operator << (ostream &outputStr, const CSampleIDWeight &sample)
{
	outputStr << sample.data.dim << "\t"; 
	for (int i=0; i<sample.data.dim; i++)
		outputStr << sample.data[i] << "\t"; 
	outputStr << sample.id << "\t"; 
	outputStr << sample.weight << "\t"; 
	outputStr << sample.calculated << endl; 

	return outputStr; 
}

bool LoadSampleFromFile(const string &file_name, vector<CSampleIDWeight> &Y)
{
	ifstream input_file(file_name.c_str(), ios::binary|ios::in);
	if (!input_file)
		return false;
	CSampleIDWeight sample; 
	read(input_file, &sample); 
	input_file.seekg(0,ios::beg); 
	input_file.seekg(0,ios::end); 
	size_t lenFile = input_file.tellg(); 
	size_t nSample = lenFile/sample.GetSize_Data(); 
	Y.resize(nSample); 
	input_file.seekg(0, ios::beg); 
	for (int i=0; i<(int)nSample; i++)	
		read(input_file, &(Y[i])); 

	input_file.close(); 
	return true; 		
}

bool SaveSampleToFile(const string &file_name, const vector<CSampleIDWeight> &Y)
{
	ofstream output_file(file_name.c_str(), ios::binary|ios::out);
	if (!output_file)
		return false; 
	for (int i=0; i<(int)Y.size(); i++)
		write(output_file, &(Y[i]) ); 
	output_file.close(); 
	return true;  
}
