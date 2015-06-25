#ifndef _MASTER_DEPLOY_HEADER_
#define _MASTER_DEPLOY_HEADER_

std::vector<CSampleIDWeight> HighestPlus1Stage_Prior( int nNode, int nInitial, CEquiEnergyModel &model);

void DispatchTuneSimulation(int nNode, int nInitial, CEquiEnergyModel &model, const CSampleIDWeight &mode, size_t simulation_length, bool save_space_flag = true, int nGroup_NSE = 0);

std::vector<CSampleIDWeight> DispatchSimulation(int nNode, int nInitial, CEquiEnergyModel &model, int simulation_length, int stage, int message_tag) ; 


#endif
