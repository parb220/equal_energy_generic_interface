#ifndef _MASTER_DEPLOY_HEADER_
#define _MASTER_DEPLOY_HEADER_

void HighestPlus1Stage(int nNode, int nInitial, CEquiEnergyModel &model); 

void DispatchHillClimbTask(int nNode, int nInitial, CEquiEnergyModel &model, int number_hill_climb);

void DispatchGMMSimulationTask(int nNode, int nInitial, CEquiEnergyModel &model, int simulation_length); 

void DispatchTuneSimulation(int nNode, int nInitial, CEquiEnergyModel &model, const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag = true);

double DispatchSimulation(int nNode, int nInitial, CEquiEnergyModel &model, size_t simulation_length, int level, int tag);


#endif
