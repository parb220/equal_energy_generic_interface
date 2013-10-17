#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>
#include <ctime>
#include <fstream>
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CEquiEnergyModel.h"
#include "CMetropolis.h"


extern "C" {
        #include "dw_math.h"
        #include "dw_switchio.h"
        #include "dw_rand.h"
}

using namespace std; 

double CEquiEnergyModel::BurnIn_RandomBlock(size_t burn_in_length)
{
	CSampleIDWeight x_new;
        int nJump =0;
        double max_posterior = current_sample.weight, bounded_log_posterior_new;
        for (int i=0; i<burn_in_length; i++)
	{
		if (metropolis->RandomBlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, 1) )
		{
			current_sample = x_new;
                        current_sample.id = (int)(time(NULL)-timer_when_started);
                        max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior;
                        nJump ++;
		}
	}
	cout << "MH Jump " << nJump << " out of " << burn_in_length << " in burning-in.\n"; 
	return max_posterior; 
}

int CEquiEnergyModel::EE_Draw_RandomBlock(size_t MH_thin)
{
        CSampleIDWeight x_new;
        int new_sample_code = NO_JUMP;
        double bounded_log_posterior_new;

        if (MakeEquiEnergyJump(x_new, current_sample))
	{
		Take_Sample_Just_Drawn_From_Storage(x_new); 
                new_sample_code = EQUI_ENERGY_JUMP;
	}
	else 
        {
                if (metropolis->RandomBlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, MH_thin))
                {
                        current_sample = x_new;
                        current_sample.id = (int)(time(NULL)-timer_when_started);
                        new_sample_code = METROPOLIS_JUMP;
                }
        }

        return new_sample_code;

}	

double CEquiEnergyModel::Simulation_Within_RandomBlock(bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new;
        int nJump =0;
        double max_posterior = current_sample.weight, bounded_log_posterior_new;
        bool if_write_file = false;
        ofstream output_file;
        if (!sample_file_name.empty() )
        {
                output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
                if (output_file)
                        if_write_file = true;
        }

        for (int i=0; i<parameter->simulation_length; i++)
        {
                for (int j=0; j<parameter->thin; j++)
		if (metropolis->RandomBlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample, parameter->THIN/parameter->thin) )
        	{
			current_sample = x_new;
               		current_sample.id = (int)(time(NULL)-timer_when_started);
               		max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior;
               		nJump ++;
        	}
		if (if_storage)
       			SaveSampleToStorage(current_sample); 
       		if (if_write_file)
      			 write(output_file, &current_sample);
       }
       if (if_write_file)
       			output_file.close();
	cout << "MH Jump " << nJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n"; 
	return max_posterior;
}

double CEquiEnergyModel::Simulation_Cross_RandomBlock(bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new;

        int nEEJump=0, nMHJump=0;
        bool if_write_file = false;
        ofstream output_file;
        if (!sample_file_name.empty() )
        {
                output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
                if (output_file)
                        if_write_file = true;
        }

        double max_posterior = current_sample.weight;
        for (int i=0; i<parameter->simulation_length; i++)
        {
                for (int j=0; j<parameter->thin; j++)
                {
                        int jump_code = EE_Draw_RandomBlock(parameter->THIN/parameter->thin);
                        if (jump_code == EQUI_ENERGY_JUMP)
                                nEEJump++;
                        else if (jump_code == METROPOLIS_JUMP)
                                nMHJump++;
                        if (jump_code == EQUI_ENERGY_JUMP || jump_code == METROPOLIS_JUMP)
                                max_posterior = max_posterior > current_sample.weight ? max_posterior : current_sample.weight;
                }

                if (if_storage)
                        SaveSampleToStorage(current_sample); 
                if (if_write_file)
                        write(output_file, &current_sample);
        }
	if (if_write_file)
                output_file.close();

        cout << "EE Jump " << nEEJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n";
        cout << "MH Jump " << nMHJump << " out of " << parameter->simulation_length *parameter->thin<< " in simulation.\n";
        return max_posterior;
}
