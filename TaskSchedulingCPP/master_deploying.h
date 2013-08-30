#ifndef _MASTER_DEPLOY_HEADER_
#define _MASTER_DEPLOY_HEADER_

double DispatchHillClimbTask(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, int number_hill_climb);

double DispatchTuneSimulation(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, size_t n_initial);

double TopDownTuningSimulation(CEquiEnergy_TState &model, const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, const CSampleIDWeight &mode, size_t period, size_t max_period);

double DispatchSimulation(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, int level, int tag);

#endif
