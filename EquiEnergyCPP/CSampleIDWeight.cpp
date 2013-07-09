#include <string>
#include "CSampleIDWeight.h"

using namespace std; 

CSampleIDWeight::CSampleIDWeight(const TDenseVector &_x, int _id, double _weight) : id(_id), weight(_weight)
{
	data.CopyContent(_x); 
}

CSampleIDWeight::CSampleIDWeight(const CSampleIDWeight &right) : id(right.id), weight(right.id)
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
	return *this; 
}

bool CSampleIDWeight::PartialCopyFrom(const CSampleIDWeight &right, unsigned int offset, size_t length)
{
	if (data.dim < offset+length || right.data.dim < offset+length)
	{
		cerr << "CSampleIDWeight::PartialCopyFrom() : index exceeding dimension" << endl; 
		return false; 
	}

	for (unsigned int i=0; i<length; i++)
		data[i+offset] = right.data[i+offset]; 	
	
	if (offset ==0 && length == data.dim && length == right.data.dim )
		weight = right. weight; 

	return true; 
}

bool CSampleIDWeight::PartialCopyFrom(unsigned int offset1, const CSampleIDWeight &right, unsigned offset2, size_t length)
{
	if (data.dim < offset1+length || right.data.dim < offset2+length )
	{
		cerr << "CSampleIDWeight::PartialCopyFrom() : index exceeding dimension" << endl; 
		return false; 
	}
	for (unsigned int i=0; i<length; i++)
		data[offset1+i] = right.data[offset2+i]; 

	if (offset1 == 0 && offset2 == 0 && length == data.dim && length == right.data.dim)
		weight = right.weight;
	
	return true; 
}

istream & read (istream & input_stream, CSampleIDWeight *x)
{
	int dim; 
	input_stream.read((char*)&(dim), sizeof(int)); 
	x->data.Resize(dim);
	double temp_x; 
	for (unsigned int i=0; i<dim; i++)
	{
		input_stream.read((char*)&temp_x, sizeof(double)); 
		x->data[i] = temp_x; 
	}
	input_stream.read((char*)&(x->id), sizeof(int)); 
	input_stream.read((char*)&(x->weight), sizeof(double)); 
	return input_stream; 
}

ostream& write(ostream & output_stream, const CSampleIDWeight *x)
{
	output_stream.write((char*)&(x->data.dim), sizeof(int)); 
	double temp_x; 
	for (unsigned int i=0; i<x->data.dim; i++)
	{
		temp_x = x->data[i]; 
		output_stream.write((char*)&temp_x, sizeof(double)); 
	}
	output_stream.write((char*)&(x->id), sizeof(int)); 
	output_stream.write((char*)&(x->weight), sizeof(double)); 
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
	return inputStr; 
}

ostream& operator << (ostream &outputStr, const CSampleIDWeight &sample)
{
	outputStr << sample.data.dim << "\t"; 
	for (int i=0; i<sample.data.dim; i++)
		outputStr << sample.data[i] << "\t"; 
	outputStr << sample.id << "\t"; 
	outputStr << sample.weight << endl; 

	return outputStr; 
}

bool LoadSampleFromFile(const string &file_name, vector<CSampleIDWeight> &Y)
{
	ifstream input_file(file_name.c_str(), ios::binary|ios::in);
	if (!input_file)
		return false;
	CSampleIDWeight sample; 
	while (!input_file.eof())
	{
		read(input_file, &sample); 
		Y.push_back(sample); 
	}	 

	input_file.close(); 
	return true; 		
}

bool SaveSampleToFile(const string &file_name, const vector<CSampleIDWeight> &Y)
{
	ofstream output_file(file_name.c_str(), ios::binary|ios::out);
	if (!output_file)
		return false; 
	for (unsigned int i=0; i<Y.size(); i++)
		write(output_file, &(Y[i]) ); 
	output_file.close(); 
	return true;  
}
