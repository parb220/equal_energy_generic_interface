#ifndef _MASTER_DEPLOY_HEADER_
#define _MASTER_DEPLOY_HEADER_

#include <mpi.h>
#include <cmath>
#include "CStorageHead.h"
#include "CEESParameter.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"
#include "mpi_parameter.h"

extern "C"
{
        #include "dw_parse_cmd.h"
}

double DispatchHillClimbTask(const vector<unsigned int> &nodePool, const CEESParameter &parameter, CStorageHead &storage, int number_hill_climb);

double DispatchTuneSimulation(const vector<unsigned int> &nodePool, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length);

double TopDownTuningSimulation(CEquiEnergyModel &model, const vector<unsigned int> &nodePool, const CEESParameter &parameter, CStorageHead &storage, const CSampleIDWeight &mode, size_t period, size_t max_period);

double DispatchSimulation(const vector<unsigned int> &nodePool, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, unsigned int level, int tag);

#endif
