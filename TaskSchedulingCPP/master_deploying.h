#ifndef _MASTER_DEPLOY_HEADER_
#define _MASTER_DEPLOY_HEADER_

void DispatchHillClimbTask(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, int number_hill_climb);

void DispatchTuneSimulation(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, size_t n_initial);


void DispatchSimulation(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, int level, int tag);


#endif
