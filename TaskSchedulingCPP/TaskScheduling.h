#ifndef _TASK_SCHEDULING_HEADER
#define _TASK_SCHEDULING_HEADER

void master_deploying(int, char **, CEquiEnergy_TState &, CEESParameter &, CStorageHead &, const CSampleIDWeight &mode );
void slave_computing(int, char **, CEquiEnergy_TState &, CEESParameter &, CStorageHead &, const CSampleIDWeight &mode);

#endif
