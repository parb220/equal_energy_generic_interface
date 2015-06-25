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
	int storage_marker;		// number of samples stored in memory for retrieve 
	string run_id; 
	int number_energy_stage;	// number of energy stages
	int number_striation; 	// number of rings
	double pee; 	// probability of equi-energy-jump
	double tN_1; 	// temperature for the highest energy stage 
	vector <double> t; 	// temperature for all stage 
	
	int highest_stage; 
	int lowest_stage; 
	int THIN;  	// thinning factor for MH samples
	int simulation_length; 
	int burn_in_length; 
public:
	CEESParameter(); 
	~CEESParameter(); 
	
	bool SetTemperature_geometric(); 
	bool SetTemperature_quadratic();  
	bool SetTemperature_polynomial(double r);

	double LogRatio_Stage(double energy_x, double energy_y, int stage) const; 
};
#endif
