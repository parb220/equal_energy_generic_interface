#ifndef _COMPUTE_NODES_
#define _COMPUTE_NODES_

void ComputeNodeLoop(CEquiEnergyModel &model);

double ExecutingTuningTask(CEquiEnergyModel &model, int stage);
double ExecutingInitialTuningTask(CEquiEnergyModel &model, int stage);
void ExecutingSimulationTask(CEquiEnergyModel &model, int stage, double scale);
void ExecutingInitialMetropolisSimulationTask(CEquiEnergyModel &model, int stage, double scale);
void ExecutingInitialIndependentSimulationTask(CEquiEnergyModel &model, int stage);

#endif
