#include "CEESParameter.h"
#include "mpi_parameter.h"

bool GetCommunicationParameter(const double *rPackage, size_t package_size, CEESParameter *parameter)
{
	parameter->simulation_length = (size_t)(rPackage[LENGTH_INDEX]); 
	parameter->burn_in_length = (size_t)(rPackage[BURN_INDEX]); 
	parameter->THIN = (size_t)(rPackage[THIN_INDEX]); 
	parameter->pee = rPackage[PEE_INDEX]; 

        // parameter->SetTemperature();
	return true; 
}

