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
	//bool SetTemperature(int nGeometricLevel, int nFinerLevel);

	double LogRatio_Level(double energy_x, double energy_y, int level) const; 

public: // parameters that are not saved, just for purpose of encapsulation
	int highest_level; 
	int lowest_level; 
	size_t thin;	// thinning factor for the entire sample
	size_t THIN;  	// thinning factor for MH samples
	size_t simulation_length; 
	size_t burn_in_length; 
	size_t max_energy_tuning_time; 
	size_t shuffle_frequency;
        size_t size_per_block;

public:
        double p_save;               // probability saving draw - actual probability of making equi-energy jump is p_save*pee
        double p_select;             // probability of selecting an independent direction to make a full size jump in that direction
	double tiny;                 // scale factor when not making a full sized jump
        int nImportanceSamples;      // number of random importance weighted previous samples to carry
	double max_energy;           // highest temperature level
	double min_ess;              // min_ess minimum desired effective sample size (default 0.80 * simulation_length?).
	int number_rings;
	int nu;                      // degrees of freedom for initial t-distribuion

	int N;                       // number to simulate for each group
	int G;                       // total number of groups is G * n_compute_cores
};
#endif
