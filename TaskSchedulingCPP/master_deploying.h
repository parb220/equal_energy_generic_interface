#ifndef _MASTER_DEPLOY_HEADER_
#define _MASTER_DEPLOY_HEADER_

void DispatchHillClimbTask(const vector<vector<int> > &nodeGroup, CEquiEnergyModel &model, int number_hill_climb);

void DispatchTuneSimulation(const vector<vector<int> > &nodeGroup, CEquiEnergyModel &model, const CSampleIDWeight &mode, size_t simulation_length);

void DispatchSimulation(const vector<vector<int> > &nodeGroup, CEquiEnergyModel &model, size_t simulation_length, int level, int tag);


#endif
