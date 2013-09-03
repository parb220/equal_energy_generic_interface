#ifndef _C_EES_PARAMETER
#define _C_EES_PARAMETER

#include <string>
#include <vector>

using namespace std; 

class CSampleIDWeight; 
bool Solve_Polynomial_Equation(vector<double> &, size_t n, double t0, double tk_1); 

class CEESParameter 
{
public:
	// The following parameter can either be determined from reading the file or command line or by default values (specified in a header file)
	string storage_dir;	// fold to store all relevant information
	size_t storage_marker;		// number of samples stored in memory for retrieve 
	string run_id; 
	size_t number_energy_level;	// number of energy levels
	double pee; 	// probability of equi-energy-jump
	double h0; 	// energy bound for the lowest energy level
	double t0; 	// temperature for the lowest energy level
	double tk_1; 	// temperature for the highest energy level

	vector <double> t; 	// temperature for all levels
public:
	CEESParameter(); 
	~CEESParameter(); 
	bool LoadParameterFromFile(string); 
	bool SaveParameterToFile(string) const;
	bool WriteSummaryFile(string) const; 
	
	bool SetEnergyBound();
	bool SetTemperature();

	double LogRatio_Level(double energy_x, double energy_y, int level) const; 

public: // parameters that are not saved, just for purpose of encapsulation
	int highest_level; 
	int lowest_level; 
	size_t deposit_frequency; 
	size_t simulation_length; 
	size_t burn_in_length; 
	size_t max_energy_tuning_time; 
	size_t shuffle_frequency;
        size_t size_per_block;
};
#endif
