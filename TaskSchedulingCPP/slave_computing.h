#ifndef _SLAVE_COMPUTING_HEADER
#define _SLAVE_COMPUTING_HEADER

bool GetCommunicationParameter(const double *, size_t, CEESParameter *); 
// CEESParameter cannot be const, because its h and t will be altered upon the received message

std::vector<int> ExecutingSimulationTask(TDenseMatrix &jump_table, CEquiEnergyModel &model, int my_rank, int group_index, int nGroup, const CSampleIDWeight &mode, int message_tag); 
// target cannot be const, because its model will be changed logically
// // storage cannot be const, because its bins will be altered constantly for deposition and drawing

#endif
