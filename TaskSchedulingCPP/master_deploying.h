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

double DispatchHillClimbTask(const vector<vector<unsigned int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, int number_hill_climb);

double DispatchTuneSimulation(const vector<vector<unsigned int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, size_t n_initial);

double TopDownTuningSimulation(CEquiEnergyModel &model, const vector<vector<unsigned int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, const CSampleIDWeight &mode, size_t period, size_t max_period);

bool ConsolidateSampleForCovarianceEstimation(const string &file_pattern, const string &variance_file); 
size_t glob(vector<string> &filename, const string &pattern); 
double DispatchSimulation(const vector<vector<unsigned int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, unsigned int level, int tag);

#endif
