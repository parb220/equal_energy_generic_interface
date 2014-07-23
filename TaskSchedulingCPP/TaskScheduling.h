#ifndef _TASK_SCHEDULING_HEADER
#define _TASK_SCHEDULING_HEADER

void master_deploying(int, int, int, CEquiEnergyModel &, const CSampleIDWeight &mode);
void slave_computing(int, int, int, CEquiEnergyModel &, const CSampleIDWeight &mode, int =10, int =10, double = 1.0);

#endif
