#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include "dw_dense_matrix.hpp"

using namespace std; 

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cerr << argv[0] << " input_filename(binary) output_filename(text)"; 
		abort(); 
	} 

	fstream iFile; 
	iFile.open(argv[1], ios::in|ios::binary); 
	if(!iFile)
	{
		cout << "Error in reading " << argv[1] << endl; 
		exit(-1); 
	}
	size_t n_blocks;
        iFile.read((char*)(&n_blocks),sizeof(size_t) );
       	vector<TDenseMatrix>blocks(n_blocks);
        for (int i=0; i<n_blocks; i++)
        {
                size_t m, n;
                iFile.read((char*)(&m),sizeof(size_t) );
                iFile.read((char*)(&n),sizeof(size_t) );
                blocks[i]=TDenseMatrix(m,n);
                for (int j=0; j<n; j++)
                	for (int k=0; k<m; k++)
                        {
                                double data;
                                iFile.read((char*)(&data),sizeof(double));
                                blocks[i](k,j)=data;
                        }
        }

        iFile.close();

	ofstream oFile; 
	oFile.open(argv[2]); 
	if (!oFile)
	{
		cout << "Error in writing " << argv[2] << endl; 
		exit(-1);
	}

	for (int i=0; i<(int)blocks.size(); i++)
	{
		oFile << blocks[i].rows << "\t" << blocks[i].cols << endl; 
		// for (int j=0; j<blocks[i].rows ; j++)
		//{
		//	for (int k=0; k<blocks[i].cols; k++)
		//		oFile << setprecision(20) << blocks[i](j,k) << "\t"; 
		//	oFile << endl; 
		//}
		oFile << setprecision(20) << blocks[i] << endl; 
	}

	oFile.close(); 
}
