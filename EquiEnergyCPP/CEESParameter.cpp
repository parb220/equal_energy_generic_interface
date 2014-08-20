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

bool CEESParameter::SetTemperature(bool if_additional)
{
        /*
 *  *     T[i] = T[i-1]+gamma^i
 *   *     gamma is determined by solving a polynomial equation 
 *    *     gamma+gamma^2+...+gamma^{K-1} = T[K-1]-T[0]; 
 *     *     */
	if (!if_additional)
		return Solve_Polynomial_Equation(t, number_energy_level+1, t0, tk_1); // t[0]=t0, t[1], ..., t[number_energy_level-1]=tk_1, t[number_energy_level] 
	else 
	{
		// Solve_Polynomial_Equation(t, number_energy_level+2, t0, tk_1); // t[0]=t0, t[1], ..., t[number_energy_level]=tk_1, t[number_energy_level+1] 
		// t.erase(t.end()-1); // remove t[number_energy_level+1], so now we have t[0]=t0, t[1], ..., t[number_energy_level]=tk_1
		t.resize(number_energy_level+1); 
		t[0] = t0; 
		t[number_energy_level] = tk_1; 
		for (int i=1; i<number_energy_level; i++) 
			t[i] = t[0] + (t[number_energy_level]-t[0])*i/(double)number_energy_level; 
	}

	/*t.resize(number_energy_level+1);  
	if (number_energy_level == 1)
	{
		t[number_energy_level] = t[0] = t0; 
		return true; 
	}
	else 
	{
		double *coefficients = new double [number_energy_level];
        	coefficients[0] = t0-tk_1;
        	for (int i=1; i<number_energy_level; i++)
                	coefficients[i]=1.0;
        	double *Z = new double [(number_energy_level-1)*2];

        	gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(number_energy_level);
        	gsl_poly_complex_solve(coefficients, number_energy_level, w, Z);

        	double gamma;
        	bool continue_flag = true;
        	for (int i=0; i<number_energy_level-1 && continue_flag; i++)
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

		t[0] = t0; 
		t[number_energy_level-1] = tk_1; 
		for (int i=1; i<number_energy_level-1; i++)
			t[i] = t[i-1]+pow(gamma, (double)i);
		t[number_energy_level] = t[number_energy_level-1]+pow(gamma, (double)number_energy_level);  
		return true; 
	}*/
}

double CEESParameter::LogRatio_Level(double original_energy_x, double original_energy_y, int level) const
{
	double log_prob_x_bounded = -original_energy_x/t[level];
        double log_prob_y_bounded = -original_energy_y/t[level];
        return log_prob_x_bounded - log_prob_y_bounded;
}

