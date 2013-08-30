#ifndef _SLAVE_COMPUTING_HEADER
#define _SLAVE_COMPUTING_HEADER

bool GetCommunicationParameter(const double *, size_t, CEESParameter &); 
// CEESParameter cannot be const, because its h and t will be altered upon the received message

bool ExecutingSimulationTask(double &min_energy, bool if_within, bool if_write_sample_file, bool if_storage, CEquiEnergy_TState &model, CStorageHead &storage, const CEESParameter &parameter, int my_rank, int group_index, size_t initialPoolSize, const CSampleIDWeight &mode, int message_tag); 
// target cannot be const, because its model will be changed logically
// // storage cannot be const, because its bins will be altered constantly for deposition and drawing
//
bool ExecutingTuningTask_BeforeSimulation(size_t period, size_t max_period, CEquiEnergy_TState &model, CStorageHead &storage, const CEESParameter &parameter, int group_index, size_t pool_size, const CSampleIDWeight &mode); 

bool ExecutingTuningTask_AfterSimulation(size_t period, size_t max_period, CEquiEnergy_TState &model, const CEESParameter &parameter, int group_index, const CSampleIDWeight &mode); 

#endif
