#include <getopt.h>
#include <fstream>
#include <iomanip>
#include "dw_dense_matrix.hpp"

using namespace std; 

int main(int argc, char **argv)
{
	static struct option long_options[] = 
	{
		{"binary file", required_argument, 0, 'b'}, 
		{"text file", required_argument, 0, 't'},
		{0, 0, 0, 0}
	};	

	int option_index = 0; 
	string binary_file, text_file; 
	while(1)
	{
		int c = getopt_long(argc, argv, "b:t:", long_options, &option_index); 
		if (c == -1)
			break; 
		switch(c)
		{
			case 'b': 
				binary_file = string(optarg); break; 
			case 't': 
				text_file = string(optarg); break; 
		}
	}

	if (binary_file.empty() || text_file.empty())
	{
		cerr << "Usage: " << argv[0] << " -b binary file -t text file.\n"; 
		abort(); 
	}

	fstream input_file; 
	input_file.open(binary_file.c_str(), ios::in | ios::binary); 
	if ( !input_file.is_open() )
	{
		cerr << "Error in opening " << input_file << " for reading.\n"; 
		abort(); 
	}
	TDenseMatrix draws; 
	draws.ReadBinary(input_file); 
	input_file.close(); 

	ofstream output_file; 
	output_file.open(text_file.c_str(), ios::out); 
	if ( !output_file.is_open() )
	{
		cerr << "Error in opening " << output_file << " for writing.\n"; 
		abort(); 
	}
	for (int i=0; i<draws.rows; i++)
	{
		for (int j=0; j<draws.cols; j++)
			output_file << setprecision(20) << draws(i,j) << "\t"; 
		output_file << endl; 
	}	
	output_file.close(); 
}
