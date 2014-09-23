#ifndef _SLAVE_COMPUTING_HEADER
#define _SLAVE_COMPUTING_HEADER

bool GetCommunicationParameter(const double *, size_t, CEESParameter *); 
// CEESParameter cannot be const, because its h and t will be altered upon the received message

void ExecutingSimulationTask(bool if_within, bool if_write_sample_file, bool if_storage, CEquiEnergyModel &model, int my_rank, int group_index, size_t initialPoolSize, const CSampleIDWeight &mode, int message_tag); 
// target cannot be const, because its model will be changed logically
// // storage cannot be const, because its bins will be altered constantly for deposition and drawing
//
bool ExecutingTuningTask_BeforeSimulation(size_t period, size_t max_period, CEquiEnergyModel &model, int group_index); 

bool ExecutingTuningTask_AfterSimulation(size_t period, size_t max_period, CEquiEnergyModel &model, int group_index); 

void ExecutingInitialTuningTask(CEquiEnergyModel &model, int node);

void ExecutingInitialSimulationTask(CEquiEnergyModel &model, int nNode, int node);

#endif
