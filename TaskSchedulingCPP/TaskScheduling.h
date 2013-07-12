#ifndef _TASK_SCHEDULING_HEADER
#define _TASK_SCHEDULING_HEADER

#include "CStorageHead.h"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"

void master_deploying(int, char **, CEquiEnergyModel &, CEESParameter &, CStorageHead &, const CSampleIDWeight &mode );
void slave_computing(int, char **, CEquiEnergyModel &, CEESParameter &, CStorageHead &, const CSampleIDWeight &mode);

#endif
