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
	lambda_1(0.0), 
	lambda(vector<double>(0)),
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
	// lambda
	if (number_energy_stage <= 0)
		return false; 
	lambda.resize(number_energy_stage+1); 
	lambda[0] = 1.0;
	if (number_energy_stage > 1)
	{
		double r = log(lambda_1)/(number_energy_stage-1); 
		for (int i=1; i<number_energy_stage+1; i++)
			lambda[i] = exp((double)i*r); 	
	}
	return true; 
}

bool CEESParameter::SetTemperature_quadratic()
{
	// lambda
	if (number_energy_stage <= 0)
		return false; 
	lambda.resize(number_energy_stage+1); 
	lambda[0] = 1.0;
	for (int i=1; i<number_energy_stage; i++)
		lambda[i] = ((double)(number_energy_stage-i) * (double)(number_energy_stage-i)) /((double)number_energy_stage*(double)number_energy_stage ); 
	lambda[number_energy_stage] = 0.1*0.1 /((double)number_energy_stage*(double)number_energy_stage );
	return true; 
}

bool CEESParameter::SetTemperature_polynomial(double r)
{
	if (number_energy_stage <= 0)
		return false; 	
	lambda.resize(number_energy_stage+1); 
	lambda[0] = 1.0; 
	double denominator = pow((double)number_energy_stage,r); 
	for (int i=1; i<number_energy_stage; i++)
		lambda[i] = pow((double)(number_energy_stage-i),r)/denominator; 
	lambda[number_energy_stage] = pow(0.1,r)/denominator; 
	return true; 
}

double CEESParameter::LogRatio_Stage(const CSampleIDWeight &x, const CSampleIDWeight &y, int stage) const 
{
	double log_prob_x_bounded = lambda[stage]*x.weight; 
	double log_prob_y_bounded = lambda[stage]*y.weight; 
        return log_prob_x_bounded - log_prob_y_bounded;
}


