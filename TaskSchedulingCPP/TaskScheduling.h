#ifndef _TASK_SCHEDULING_HEADER
#define _TASK_SCHEDULING_HEADER

void master_deploying(int, CEquiEnergyModel &);
void DispatchTuneSimulation(int, CEquiEnergyModel &, size_t, bool);
void slave_computing(CEquiEnergyModel &);

#endif
