#include "slave_computing.h"

using namespace std;

void slave_computing(int argc, char **argv, CEquiEnergyModel &model, CEESParameter &parameter, CStorageHead &storage, const CSampleIDWeight &mode) 
{
	int my_rank, nCPU; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &nCPU);
	
	MPI_Status status; 
	
	double *rPackage = new double [N_MESSAGE]; 
	double *sPackage = new double [N_MESSAGE];    
	
	double min_energy; 
	bool if_within, if_write_file, if_storage;   

	while (1)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
		if (status.MPI_TAG == END_TAG)
		{
			delete []rPackage; 
			delete []sPackage; 
			exit(0); 
		}		
		else if (status.MPI_TAG == HILL_CLIMB_TAG)
			model.HillClimb_CSMINWEL(rPackage[LENGTH_INDEX], storage, parameter); 
			// model.HillClimb_NPSOL(rPackage[LENGTH_INDEX], storage, parameter); 
		else if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST || status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND || status.MPI_TAG == SIMULATION_TAG) 
		{	
			model.energy_level = (unsigned int)(rPackage[LEVEL_INDEX]);
			if (!GetCommunicationParameter(rPackage, N_MESSAGE, parameter))
			{
				cout << "GetCommunicationParameter() : Error occurred.\n"; 
				abort(); 
			}
			model.h_bound = parameter.h[model.energy_level]; 
			model.t_bound = parameter.t[model.energy_level];

			if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST)
				if_within = true; 
			else 
				if_within = false; 
	
			if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST || status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND)
				if_write_file = true; 
			else 
				if_write_file = false; 

			if (status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND || status.MPI_TAG == SIMULATION_TAG)
				if_storage = true;
			else 
				if_storage = false; 	

			bool simulation_flag = ExecutingSimulationTask(min_energy, if_within, if_write_file, if_storage, model, storage, parameter, my_rank, 2*(size_t)nCPU, mode, status.MPI_TAG); 

			if (!simulation_flag)
			{
				cerr << "Error in simulation.\n"; 
				abort(); 
			}
			sPackage[LEVEL_INDEX] = model.energy_level;
			sPackage[H0_INDEX] = min_energy; 	
			MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
		}
	}
}

bool GetCommunicationParameter(const double *rPackage, size_t package_size, CEESParameter &parameter)
{
	parameter.simulation_length = (size_t)(rPackage[LENGTH_INDEX]); 
	parameter.burn_in_length = (size_t)(rPackage[BURN_INDEX]); 
	parameter.deposit_frequency = (size_t)(rPackage[FREQ_INDEX]); 

        parameter.h0 = rPackage[H0_INDEX];
        parameter.SetEnergyBound();
        parameter.SetTemperature();
	return true; 
}

