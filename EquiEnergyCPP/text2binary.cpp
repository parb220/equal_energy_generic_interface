#include <fstream>
#include <getopt.h>
#include "CSampleIDWeight.h"

using namespace std; 

int main(int argc, char **argv)
{
	static struct option long_options[] = 
	{
		{"text file", required_argument, 0, 't'},
		{"binary file", required_argument, 0, 'b'}, 
		{0, 0, 0, 0}
	}; 
	
	int option_index=0; 
	string text_file, binary_file; 

	while (1)
        {
                int c = getopt_long(argc, argv, "t:b:", long_options, &option_index);
                if (c == -1)
                        break;
		switch(c)
		{
			case 't':
				text_file = string(optarg); break; 
			case 'b': 
				binary_file = string(optarg); break; 
			default: 
				break; 
		}
	}
	if (text_file.empty() || binary_file.empty())
	{
		cerr << "Usage: " << argv[0] << "-t text file -b binary file.\n"; 
		abort(); 
	}

	CSampleIDWeight sample; 
	ifstream input; 
	input.open(text_file.c_str(), ios::in); 
	if (!input)
	{
		cerr << "Error opening " << text_file << endl; 
		abort(); 
	}
	input >> sample; 
	input.close(); 

	ofstream output(binary_file.c_str(), ios::out|ios::binary); 
	if (!output)
	{
		cerr << "Error opening " << binary_file << endl; 
		abort(); 
	}
	write(output, &(sample)); 
	output.close(); 
}
