#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include "CSampleIDWeight.h"

using namespace std; 

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cerr << argv[0] << " input_filename(binary) output_filename(text)"; 
		abort(); 
	} 
	CSampleIDWeight one_sample; 
	vector <CSampleIDWeight> sample;
	if (!LoadSampleFromFile(string(argv[1]), sample) )
	{
		cerr << "Error in LoadSampleFromFile() : " << argv[1] << endl; 
		abort(); 
	}

	ofstream oFile; 
	oFile.open(argv[2]); 
	if (!oFile)
	{
		cout << "Error in writing " << argv[2] << endl; 
		exit(-1);
	}

	for (unsigned int i=0; i<sample.size(); i++)
		oFile << setprecision(20) << sample[i]; 

	oFile.close(); 
}
