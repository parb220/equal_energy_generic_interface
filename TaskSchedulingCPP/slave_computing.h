#ifndef _SLAVE_COMPUTING_HEADER
#define _SLAVE_COMPUTING_HEADER

void slave_computing(CEquiEnergyModel &model);

void ExecutingInitialIndependentSimulationTask(CEquiEnergyModel &model, int node);

double ExecutingInitialTuningTask(CEquiEnergyModel &model);

void ExecutingInitialMetropolisSimulationTask(CEquiEnergyModel &model,double scale, int node);

double ExecutingTuningTask(CEquiEnergyModel &model, int level, double Kplus);

void ExecutingSimulationTask(CEquiEnergyModel &model, int node, double scale);

#endif
