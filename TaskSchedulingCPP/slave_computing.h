#include <mpi.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include "CEquiEnergyModel.h"
#include "CStorageHead.h"
#include "CEESParameter.h"
#include "CMetropolis.h"
#include "mpi_parameter.h"

bool GetCommunicationParameter(const double *, size_t, CEESParameter &); 
// CEESParameter cannot be const, because its h and t will be altered upon the received message

bool ExecutingSimulationTask(double &min_energy, bool if_within, bool if_write_sample_file, bool if_storage, CEquiEnergyModel &model, CStorageHead &storage, const CEESParameter &parameter, unsigned int my_rank, unsigned int group_index, size_t initialPoolSize, const CSampleIDWeight &mode, int message_tag); 
// target cannot be const, because its model will be changed logically
// // storage cannot be const, because its bins will be altered constantly for deposition and drawing
//
bool ExecutingTuningTask_BeforeSimulation(size_t period, size_t max_period, CEquiEnergyModel &model, CStorageHead &storage, const CEESParameter &parameter, unsigned int group_index, size_t pool_size, const CSampleIDWeight &mode); 

bool ExecutingTuningTask_AfterSimulation(size_t period, size_t max_period, CEquiEnergyModel &model, const CEESParameter &parameter, unsigned int group_index, const CSampleIDWeight &mode); 
