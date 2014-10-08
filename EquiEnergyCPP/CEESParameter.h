#ifndef _C_EES_PARAMETER
#define _C_EES_PARAMETER

#include <string>
#include <vector>

using namespace std; 

class CSampleIDWeight; 
//bool Solve_Polynomial_Equation(vector<double> &, size_t n, double t0, double tk_1); 

class CEESParameter 
{
public:
	// The following parameter can either be determined from reading the file or command line or by default values (specified in a header file)
	string storage_dir;	// fold to store all relevant information
	size_t storage_marker;		// number of samples stored in memory for retrieve 
	string run_id; 
	size_t number_energy_level;	// number of energy levels
	double pee; 	// probability of equi-energy-jump

public:
	CEESParameter(); 
	~CEESParameter(); 
	bool LoadParameterFromFile(string); 
	bool SaveParameterToFile(string) const;
	bool WriteSummaryFile(string) const; 

public: // parameters that are not saved, just for purpose of encapsulation
	int highest_level; 
	int lowest_level; 
	size_t simulation_length; 


public:
        double expected_block_size;  // p_select = expected_block_size/NParameters.  If less than or equal to zero, p_select = 1.0.
        double p_select;             // probability of selecting an independent direction to make a full size jump in that direction
        double p_save;               // probability saving draw.  must be greater than pee.
	double pee_divided_psave;    // pee/p_save
	double tiny;                 // scale factor when not making a full sized jump
        int nImportanceSamples;      // number of random importance weighted previous samples to carry
	double max_energy;           // highest temperature level
	double min_ess;              // min_ess minimum desired effective sample size
	int number_rings;            // number of rings at each level
	int nu;                      // degrees of freedom for initial t-distribuion
	double geometric;            // if positive, multiplier for temperatures, if one or greater, is set by number of temperature levels

	int desired_G;               // desired number of groups, actual number of groups will be at least his
	int N;                       // number to simulate for each group
	int Gn;                      // total number of groups is G * n_compute_cores
};
#endif
