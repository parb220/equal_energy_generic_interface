#include <fstream>
#include <cmath>

#include "prcsn.h"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"

using namespace std; 

CEESParameter::CEESParameter() : 
	storage_dir(string()), 
	storage_marker(0), 
	run_id(string()), 
	number_energy_stage(0),
	number_striation(0), 
	pee(0.0),
	tN_1(0.0), 
	t(vector<double>(0)),
	highest_stage(0), 
	lowest_stage(0), 
	THIN(0), 
	simulation_length(0), 
	burn_in_length(0)
{}

CEESParameter::~CEESParameter()
{}		

bool CEESParameter::SetTemperature_geometric()
// t[i+1] = t[i] * r
{
	if (number_energy_stage < 0)
		return false; 

	t.resize(number_energy_stage+1); 
	t[0] = 1.0;  
	if (number_energy_stage > 1) 
	{
		double r=exp(log(tN_1)/double(number_energy_stage-1)); 
		for (int i=1; i<(int)(t.size()); i++)
			t[i] = t[i-1] * r; 
	}
	else
		t[1] = tN_1; 
	return true; 
}

bool CEESParameter::SetTemperature_quadratic()
{
	if (number_energy_stage < 0)
		return false; 

	t.resize(number_energy_stage+1);
	for (int i=0; i<number_energy_stage; i++)
		t[i] = (double)number_energy_stage*(double)number_energy_stage/((double)(number_energy_stage-i) * (double)(number_energy_stage-i)); 

	t[number_energy_stage] = PLUS_INFINITY; 
	return true; 
}

bool CEESParameter::SetTemperature_polynomial(double r)
{
	if (number_energy_stage < 0)
		return false; 
	t.resize(number_energy_stage+1);
	for (int i=0; i<number_energy_stage; i++)
		t[i] = pow((double)number_energy_stage,r)/pow((double)(number_energy_stage-i),r); 

	t[number_energy_stage] = PLUS_INFINITY; 
}

double CEESParameter::LogRatio_Stage(double original_energy_x, double original_energy_y, int stage) const
{
	double log_prob_x_bounded = -original_energy_x/t[stage];
        double log_prob_y_bounded = -original_energy_y/t[stage];
        return log_prob_x_bounded - log_prob_y_bounded;
}

