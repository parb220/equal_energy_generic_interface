#include <fstream>
#include <getopt.h>
#include <vector>
#include "CSampleIDWeight.h"

using namespace std; 

int main(int argc, char **argv)
{
	static struct option long_options[] = 
	{
		{"number of samples", required_argument, 0, 'n'}, 
		{"dimension of samples", required_argument, 0, 'd'},
		{"text file", required_argument, 0, 't'},
		{"binary file", required_argument, 0, 'b'}, 
		{0, 0, 0, 0}
	}; 
	
	int option_index=0; 
	string text_file, binary_file; 
	int nDraws=-1, dim=-1; 

	while (1)
        {
                int c = getopt_long(argc, argv, "t:b:n:d:", long_options, &option_index);
                if (c == -1)
                        break;
		switch(c)
		{
			case 't':
				text_file = string(optarg); break; 
			case 'b': 
				binary_file = string(optarg); break; 
			case 'n':
				nDraws = atoi(optarg); break; 
			case 'd': 
				dim = atoi(optarg); break; 
			default: 
				break; 
		}
	}
	
	if (text_file.empty() || binary_file.empty() || nDraws < 0 || dim < 0)
	{
		cerr << "Usage: " << argv[0] << "-t text file -b binary file -n number of samples -d dimension of samples.\n"; 
		abort(); 
	}

 	vector<CSampleIDWeight> sample(nDraws); 
	TDenseVector data(dim+2); 
	ifstream input; 
	input.open(text_file.c_str(), ios::in); 
	if (!input)
	{
		cerr << "Error opening " << text_file << endl; 
		abort(); 
	}
	for (int i=0; i<nDraws; i++)
	{
		input >> data; 
		sample[i].data.CopyContent(data.SubVector(2,dim+1)); 
		sample[i].weight = data[0]; 
	}
	input.close(); 

	ofstream output(binary_file.c_str(), ios::out|ios::binary); 
	if (!output)
	{
		cerr << "Error opening " << binary_file << endl; 
		abort(); 
	}
	for (int i=0; i<nDraws; i++)
		write(output, &(sample[i])); 
	output.close(); 
}
