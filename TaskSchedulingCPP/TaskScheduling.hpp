#ifndef TASK_SCHEDULING_HEADER
#define TASK_SCHEDULING_HEADER

void master_deploying(const int, int, int, CEquiEnergyModel &, const CSampleIDWeight &mode, int nGroup_NSE=0);
void slave_computing(const int, CEquiEnergyModel &, const CSampleIDWeight &mode);
#endif
