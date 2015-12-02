#ifndef MASTER_DEPLOY_HEADER_
#define MASTER_DEPLOY_HEADER_

#include "option.hpp"

std::vector<CSampleIDWeight> HighestPlus1Stage_Prior(double *, double *, const int, int nNode, int nInitial, CEquiEnergyModel &model, std::ofstream &);

void DispatchTuneSimulation(double *, double *, const int, int nNode, int nInitial, CEquiEnergyModel &model, const CSampleIDWeight &mode, int simulation_length, const Diagnosis &option=OPT_ESS, bool save_space_flag = true);

std::vector<CSampleIDWeight> DispatchSimulation(double *, double *, const int, int nNode, int nInitial, CEquiEnergyModel &model, int simulation_length, int stage, int message_tag, std::ofstream &) ; 


#endif
