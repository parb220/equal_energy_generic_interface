#ifndef _MASTER_DEPLOY_HEADER_
#define _MASTER_DEPLOY_HEADER_

void DispatchHillClimbTask(int nNode, CEquiEnergyModel &model, int number_hill_climb);

void DispatchGMMSimulationTask(int nNode, CEquiEnergyModel &model, int simulation_length); 

void DispatchTuneSimulation(int, CEquiEnergyModel &, size_t, bool);

double DispatchSimulation(int nNode, int nInitial, CEquiEnergyModel &model, size_t simulation_length, int level, int tag);

void DispatchInitialSimulation(int nNode, CEquiEnergyModel &model);

#endif
