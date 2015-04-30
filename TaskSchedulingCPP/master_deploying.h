#ifndef _MASTER_DEPLOY_HEADER_
#define _MASTER_DEPLOY_HEADER_

std::vector<CSampleIDWeight> HighestPlus1Stage(int nNode, int nInitial, CEquiEnergyModel &model); 
std::vector<CSampleIDWeight> HighestPlus1Stage_Prior( int nNode, int nInitial, CEquiEnergyModel &model);

void DispatchHillClimbTask(int nNode, int nInitial, CEquiEnergyModel &model, int number_hill_climb);

void DispatchGMMSimulationTask(int nNode, int nInitial, CEquiEnergyModel &model, int simulation_length); 

void DispatchTuneSimulation(int nNode, int nInitial, CEquiEnergyModel &model, const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag = true);

std::vector<CSampleIDWeight> DispatchSimulation(int nNode, int nInitial, CEquiEnergyModel &model, int simulation_length, int stage, int message_tag); 


#endif
