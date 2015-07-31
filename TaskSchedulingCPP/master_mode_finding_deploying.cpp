#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <mpi.h>
#include "CEquiEnergyModel.h"
#include "mpi_parameter.h"
#include "storage_parameter.h"

using namespace std;

void master_mode_finding_deploying(int nNode, int nInitial, CEquiEnergyModel &model, const CSampleIDWeight &mode, double pctl)
{	
	double *sPackage= new double [N_MESSAGE];
	double *rPackage= new double [N_MESSAGE];   
	MPI_Status status;

	// binning 
        model.storage->binning_equal_size(0, model.parameter->number_striation, 1.0); 

	// Starting points
	vector<CSampleIDWeight> samples = model.storage->DrawAllSample(0), start_points; 
	if (!samples.empty())
		start_points = model.Initialize_MostDistant_WithinPercentile(samples, nInitial, 0, pctl); 
	if (samples.empty() || start_points.empty())
	{
		start_points.resize(nInitial); 
		for (int i=0; i<nInitial; i++)
			start_points[i]  = mode; 
	}

	stringstream convert;
        ofstream output_file;
        convert.str(string());
        convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << "0";
        string start_point_file = model.parameter->storage_dir + convert.str();
        output_file.open(start_point_file.c_str(), ios::binary|ios::out);
        if (!output_file)
        {
        	cerr << "Error in writing to " << start_point_file << endl;
        	abort();
        }
	for (int i=0; i<nInitial; i++)
       		write(output_file, &(start_points[i]));
       	output_file.close();

	// send out hill_climb jobs
	vector<int> available_node(nNode-1);
        for (int i=0; i<(int)(available_node.size()); i++)
                available_node[i] = i+1;

	int iInitial =0; 
        while (iInitial < nInitial)
        {
               	sPackage[GROUP_INDEX] = iInitial;
		sPackage[LENGTH_INDEX] = 1; 
                if (available_node.empty())
                {
                	MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);
                       	available_node.push_back(status.MPI_SOURCE);
                }
                MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, available_node.back(), HILL_CLIMB_TAG, MPI_COMM_WORLD);
                available_node.pop_back(); 
                iInitial ++;
        }

        for (int j=0; j<(nNode-1)-(int)available_node.size(); j++)
                MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status);

	// tell all the slaves to exit by sending an empty messag with 0 simulation length 
	for (int i=1; i<nNode; i++)
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
        delete [] sPackage;
        delete [] rPackage;
}

