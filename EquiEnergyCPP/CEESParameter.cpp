#include <fstream>
#include <gsl/gsl_poly.h>
#include <cmath>

#include "CSampleIDWeight.h"
#include "CEESParameter.h"

using namespace std; 

CEESParameter::CEESParameter() : 
	storage_dir(string()), 
	storage_marker(0), 
	run_id(0), 
	number_energy_level(0), 
	pee(0.0),
	h0(0.0), 
	hk_1(0.0),
	t0(0.0), 
	c_factor(0.0), 
	h(vector<double>(0)), 
	t(vector<double>(0)),
	highest_level(0), 
	lowest_level(0), 
	deposit_frequency(0), 
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
	oFile.write((char*)(&run_id), sizeof(int)); 
	oFile.write((char*)(&number_energy_level), sizeof(size_t)); 
	oFile.write((char*)(&pee), sizeof(double)); 
	oFile.write((char*)(&h0), sizeof(double)); 
	oFile.write((char*)(&hk_1), sizeof(double)); 
	oFile.write((char*)(&t0), sizeof(double)); 
	oFile.write((char*)(&c_factor), sizeof(double));
	
	for (unsigned int i=0; i<=number_energy_level; i++)
		oFile.write((char*)(&(h[i])), sizeof(double));

	for (unsigned int i=0; i<=number_energy_level; i++)
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
	char *storage_dir_array = new char [name_len+1]; 
	iFile.read(storage_dir_array, sizeof(char)*name_len); 
	storage_dir = string(storage_dir_array); 
	delete [] storage_dir_array; 
        
	iFile.read((char*)(&storage_marker), sizeof(size_t));
        iFile.read((char*)(&run_id), sizeof(int));
        iFile.read((char*)(&number_energy_level), sizeof(size_t));
        iFile.read((char*)(&pee), sizeof(double));
        iFile.read((char*)(&h0), sizeof(double));
        iFile.read((char*)(&hk_1), sizeof(double));
        iFile.read((char*)(&t0), sizeof(double));
        iFile.read((char*)(&c_factor), sizeof(double));
	
	h.resize(number_energy_level+1);
	for (unsigned int i=0; i<=number_energy_level; i++) 
		iFile.read((char*)(&(h[i])), sizeof(double));
        
	t.resize(number_energy_level+1);
	for (unsigned int i=0; i<=number_energy_level; i++) 
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
	oFile << "Energy Thresholds:"; 
	for (unsigned int i=0; i<=number_energy_level; i++)
		oFile << "\t" << h[i]; 
	oFile << endl; 
	oFile << "Temperatures:"; 
	for (unsigned int i=0; i<=number_energy_level; i++)
		oFile << "\t" << t[i]; 
	oFile << endl; 
	oFile << "Prob Equi-Jump:\t" << pee << endl;  
	oFile.close(); 
	return true; 
}

bool CEESParameter::SetEnergyBound()
{
        /*
 *  *     H[i] = H[i-1]+gamma^i
 *   *     gamma is determined by solving a polynomial equation 
 *    *     gamma+gamma^2+...+gamma^{K-1} = H[K-1]-H[0]; 
 *     *     */
	h.resize(number_energy_level+1);  
	if (number_energy_level == 1)
	{
		h[number_energy_level] = h[0] = h0; 
		return true; 
	}
	else 
	{
		double *coefficients = new double [number_energy_level];
        	coefficients[0] = h0-hk_1;
        	for (unsigned int i=1; i<number_energy_level; i++)
                	coefficients[i]=1.0;
        	double *Z = new double [(number_energy_level-1)*2];

        	gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(number_energy_level);
        	gsl_poly_complex_solve(coefficients, number_energy_level, w, Z);

        	double gamma;
        	bool continue_flag = true;
        	for (unsigned int i=0; i<number_energy_level-1 && continue_flag; i++)
        	{
                	if (Z[2*i]>0 && abs(Z[2*i+1]) <= 1.0e-6)
                	{
                        	gamma = Z[2*i];
                        	continue_flag = false;
                	}
        	}
        	delete [] Z;
        	delete [] coefficients;
        	if (continue_flag)
                	return false;

		h[0] = h0; 
		h[number_energy_level-1] = hk_1; 
		for (unsigned int i=1; i<number_energy_level-1; i++)
			h[i] = h[i-1]+pow(gamma, (double)i);
		h[number_energy_level] = h[number_energy_level-1]+pow(gamma, (double)number_energy_level);  
		return true; 
	}
}

bool CEESParameter::SetTemperature()
{
	t.resize(h.size()); 
	t[0] = t0; 
	for (unsigned int i=1; i<=number_energy_level; i++)
		t[i] = t[i-1]+(h[i]-h[i-1])/c_factor; 
	return true;
}

unsigned int CEESParameter::EnergyIndex(double energy) const
//  h[0] < h[1] < h[2] < ... < h[k-1] < h[k] < infinity
//  h[0] should be the lower bound of all possible energy values; 
//  that is h[0] <= energy < infinity
{
	for (unsigned int j=1; j<number_energy_level; j++)
	{
		if (energy < h[j])
			return j-1;
	}
	return number_energy_level-1; 
}

double CEESParameter::LogRatio_Level(double original_energy_x, double original_energy_y, unsigned int level) const
{
	double log_prob_x_bounded = -(original_energy_x > h[level] ? original_energy_x : h[level])/t[level];
        double log_prob_y_bounded = -(original_energy_y > h[level] ? original_energy_y : h[level])/t[level];
        return log_prob_x_bounded - log_prob_y_bounded;
}

