#include <cmath>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "dw_rand.h"
#include "dw_matrix.h"
#include "CSampleIDWeight.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"
#include "storage_parameter.h"
#include "mpi_parameter.h"
#include "master_deploying.h"

using namespace std; 

vector<string> glob(const string &pattern); 

std::vector<CSampleIDWeight> HighestPlus1Stage_Prior(int nNode, int nInitial, CEquiEnergyModel &model)
{
	return DispatchSimulation(nNode, nInitial, model, model.parameter->simulation_length,  model.parameter->number_energy_stage, SIMULATION_PRIOR_TAG);
}
