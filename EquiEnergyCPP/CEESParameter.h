#ifndef _C_EES_PARAMETER
#define _C_EES_PARAMETER

#include <string>
#include <vector>

using namespace std; 

class CSampleIDWeight; 

class CEESParameter 
{
public:
	// The following parameter can either be determined from reading the file or command line or by default values (specified in a header file)
	string storage_dir;	// fold to store all relevant information
	size_t storage_marker;		// number of samples stored in memory for retrieve 
	string run_id; 
	size_t number_energy_stage;	// number of energy stages
	size_t number_striation; 	// number of rings
	double pee; 	// probability of equi-energy-jump
	double t0; 	// temperature for the lowest energy stage 
	double tN_1; 	// temperature for the highest energy stage 

	vector <double> t; 	// temperature for all stage 
public:
	CEESParameter(); 
	~CEESParameter(); 
	bool LoadParameterFromFile(string); 
	bool SaveParameterToFile(string) const;
	bool WriteSummaryFile(string) const; 
	
	bool SetEnergyBound();
	bool SetTemperature(); 

	double LogRatio_Stage(double energy_x, double energy_y, int stage) const; 

public: // parameters that are not saved, just for purpose of encapsulation
	int highest_stage; 
	int lowest_stage; 
	size_t thin;	// thinning factor for the entire sample
	size_t THIN;  	// thinning factor for MH samples
	size_t simulation_length; 
	size_t burn_in_length; 
	size_t max_energy_tuning_time; 
	size_t shuffle_frequency;
        size_t size_per_block;
};
#endif
