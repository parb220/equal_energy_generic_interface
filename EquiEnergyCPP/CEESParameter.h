#ifndef _C_EES_PARAMETER
#define _C_EES_PARAMETER

#include <string>
#include <vector>
#include "CSampleIDWeight.h"

using namespace std; 

class CEESParameter 
{
public:
	// The following parameter can either be determined from reading the file or command line or by default values (specified in a header file)
	string storage_dir;	// fold to store all relevant information
	size_t storage_marker;		// number of samples stored in memory for retrieve 
	int run_id; 
	size_t number_energy_level;	// number of energy levels
	double pee; 	// probability of equi-energy-jump
	double h0; 	// energy bound for the lowest energy level
	double hk_1; 	// energy bound for the highest energy level
	double t0; 	// temperature for the lowest energy level
	double c_factor;	// C factor to determine temperatures based on energy bounds, (H[i]-H[i-1])/(T[i]-T[i-1])=C

	vector <double> h, t; 	// energy-bound and temperature for all levels
public:
	CEESParameter(); 
	~CEESParameter(); 
	bool LoadParameterFromFile(string); 
	bool SaveParameterToFile(string) const;
	bool WriteSummaryFile(string) const; 
	
	bool SetEnergyBound();
	bool SetTemperature();

	unsigned int BinIndex_Start(unsigned int level) const { return level*number_energy_level; }
	unsigned int BinIndex_End(unsigned int level) const { return (level+1)*number_energy_level-1; }
	unsigned int EnergyIndex(double energy) const; 
	unsigned int BinIndex(double energy, unsigned int level) const { return EnergyIndex(energy)+BinIndex_Start(level); }
	double LogRatio_Level(double energy_x, double energy_y, unsigned int level) const; 

public: // parameters that are not saved, just for purpose of encapsulation
	unsigned int highest_level; 
	unsigned int lowest_level; 
	size_t deposit_frequency; 
	size_t simulation_length; 
	size_t burn_in_length; 
	size_t max_energy_tuning_time; 
	size_t shuffle_frequency;
        size_t size_per_block;
};
#endif
