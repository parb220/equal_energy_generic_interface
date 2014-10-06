#include <fstream>
#include <cmath>

#include "CSampleIDWeight.h"
#include "CEESParameter.h"

using namespace std; 

CEESParameter::CEESParameter() : 
	storage_dir(string()), 
	storage_marker(0), 
	run_id(string()), 
	number_energy_level(0), 
	pee(0.0),
	t0(0.0), 
	tk_1(0.0), 
	t(vector<double>(0)),
	highest_level(0), 
	lowest_level(0), 
	thin(0),
	THIN(0), 
	simulation_length(0), 
	burn_in_length(0),
	max_energy_tuning_time(0),
	shuffle_frequency(0),
	size_per_block(0)
{}

CEESParameter::~CEESParameter()
{}		

bool CEESParameter::SaveParameterToFile(string file_name) const
{
	fstream oFile(file_name.c_str(), ios::out | ios::binary); 
	if (!oFile)
		return false;
	size_t storage_dir_name_len = storage_dir.length();  
	oFile.write((char*)(&storage_dir_name_len), sizeof(size_t)); 
	oFile.write((char*)(storage_dir.c_str()), sizeof(char)*storage_dir_name_len); 

	oFile.write((char*)(&storage_marker), sizeof(size_t)); 
	size_t run_id_len = run_id.length(); 
	oFile.write((char*)(&run_id_len), sizeof(size_t)); 
	oFile.write((char*)(run_id.c_str()), sizeof(char)*run_id_len); 
	oFile.write((char*)(&number_energy_level), sizeof(size_t)); 
	oFile.write((char*)(&pee), sizeof(double)); 
	oFile.write((char*)(&t0), sizeof(double)); 
	oFile.write((char*)(&tk_1), sizeof(double)); 
	
	for (int i=0; i<=number_energy_level; i++)
	 	oFile.write((char*)(&(t[i])), sizeof(double)); 

	oFile.close(); 
	return true; 
}

bool CEESParameter::LoadParameterFromFile(string file_name)
{
        fstream iFile(file_name.c_str(), ios::in | ios::binary);
        if (!iFile)
                return false;

	size_t name_len; 
	iFile.read((char*)(&name_len), sizeof(size_t)); 
	char *storage_dir_array = new char[name_len+1]; 
	iFile.read(storage_dir_array, sizeof(char)*name_len); 
	storage_dir = string(storage_dir_array); 
	delete [] storage_dir_array; 
        
	iFile.read((char*)(&storage_marker), sizeof(size_t));
	size_t run_id_len; 
	iFile.read((char*)(&run_id_len), sizeof(size_t)); 
	char *run_id_array = new char[run_id_len+1]; 
	iFile.read(run_id_array, sizeof(char)*run_id_len); 
	run_id = string(run_id_array); 
        delete [] run_id_array; 

	iFile.read((char*)(&number_energy_level), sizeof(size_t));
        iFile.read((char*)(&pee), sizeof(double));
        iFile.read((char*)(&t0), sizeof(double));
	iFile.read((char*)(&tk_1), sizeof(double)); 
	
	t.resize(number_energy_level+1);
	for (int i=0; i<=number_energy_level; i++) 
		iFile.read((char*)(&(t[i])), sizeof(double));

	iFile.close(); 
	return true; 
}

bool CEESParameter::WriteSummaryFile(string file_name) const
{
	ofstream oFile; 
	oFile.open(file_name.c_str(), ios::out); 
	if (!oFile) 
		return false; 
	oFile << "Storage Marker:\t" << storage_marker << endl; 
	oFile << "Number of Energy Levels:\t" << number_energy_level << endl; 
	oFile << "Temperatures:"; 
	for (int i=0; i<=number_energy_level; i++)
		oFile << "\t" << t[i]; 
	oFile << endl; 
	oFile << "Prob Equi-Jump:\t" << pee << endl;  
	oFile.close(); 
	return true; 
}

bool CEESParameter::SetTemperature(int nLevel)
// t[i+1] = t[i] * r
{
	number_energy_level = nLevel; 
	t.resize(number_energy_level+1); 
	double r=exp(log(tk_1/t0)/double(nLevel-1)); 
	t[0] = t0; 
	for (int i=1; i<(int)(t.size()); i++)
		t[i] = t[i-1] * r; 
	return true; 
}

bool CEESParameter::SetTemperature(int nGeometricLevel, int nFinerLevel)
// (t[i+1]-t[i])/(t[i]-t[i-1]) = r
{
	if (nGeometricLevel) // nGeometricLevel != 0
	{
		number_energy_level = nGeometricLevel + nFinerLevel; 
		t.resize(number_energy_level+1);
	
		vector<double> geometric_t(nGeometricLevel+1, 0.0); 
		Solve_Polynomial_Equation(geometric_t, nGeometricLevel+1, t0, tk_1); 

		vector<double> finer_t(nFinerLevel+3, 0.0); 
		Solve_Polynomial_Equation(finer_t, nFinerLevel+3, geometric_t[0], geometric_t[1]); 
		for (int i=0; i<nFinerLevel+1; i++)
			t[i] = finer_t[i]; 
		for (int i=0; i<nGeometricLevel; i++)
			t[nFinerLevel+i+1] = geometric_t[i+1]; 
		
		/*vector<double> finer_t(nFinerLevel, 0.0); 
		for (int i=0; i<nFinerLevel; i++)
			finer_t[i] = geometric_t[0] + (geometric_t[1]-geometric_t[0])*(double)(i+1.0)/(double)(nFinerLevel+1.0); 

		t[0] = geometric_t[0]; 
		for (int i=0; i<nFinerLevel; i++)
			t[i+1] = finer_t[i]; 
		for (int i=1; i<=nGeometricLevel; i++)
			t[nFinerLevel+i] = geometric_t[i];*/ 
	}
	else if ( nGeometricLevel == 0 && nFinerLevel != 0)
	{
		number_energy_level = nFinerLevel+2; 
		t.resize(number_energy_level+1); 
		Solve_Polynomial_Equation(t, number_energy_level+1, t0, tk_1); 
		
		/*t[0] = t0; 
		t[nFinerLevel-1] = tk_1; 
		for (int i=1; i<nFinerLevel-1; i++)
			t[i] = t[0] + (t[nFinerLevel-1]-t[0])*(double)i/(double)(nFinerLevel-1.0); 
	 	t[nFinerLevel] = t[nFinerLevel-1] + (t[nFinerLevel-1]-t[0])/(double)(nFinerLevel-1.0); */
	} 
	else // nGeometricLevel == 0 && nFinerLevel == 0
	{
		number_energy_level = 1;
		t.resize(1); 
		t[0]  = t0;  
	}
	return true; 
}

double CEESParameter::LogRatio_Level(double original_energy_x, double original_energy_y, int level) const
{
	double log_prob_x_bounded = -original_energy_x/t[level];
        double log_prob_y_bounded = -original_energy_y/t[level];
        return log_prob_x_bounded - log_prob_y_bounded;
}

