#ifndef _TASK_SCHEDULING_HEADER
#define _TASK_SCHEDULING_HEADER

void master_deploying(int, int, CEquiEnergyModel &, const CSampleIDWeight &mode, int nGroup_NSE=0);
void slave_computing(CEquiEnergyModel &, const CSampleIDWeight &mode);
#endif
