#ifndef _TASK_SCHEDULING_HEADER
#define _TASK_SCHEDULING_HEADER

void master_deploying(int, int, CEquiEnergyModel &, const CSampleIDWeight &mode);
void slave_computing(int, int, int, CEquiEnergyModel &, const CSampleIDWeight &mode, int =10, int =10, double = 1.0);

void slave_mode_finding_computing(int n_initial, CEquiEnergyModel &model, const CSampleIDWeight &mode, int optimizationN, int perturbationN, double perturbationS); 
void master_mode_finding_deploying(int nNode, int nInitial, CEquiEnergyModel &model, const CSampleIDWeight &mode, double prctl=0.10); 

void DispatchSimulation_PriorProbability(int nNode, int simulation_length); 

#endif
