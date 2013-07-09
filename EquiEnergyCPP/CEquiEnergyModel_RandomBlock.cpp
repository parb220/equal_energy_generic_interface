#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>
#include <ctime>
#include "CEquiEnergyModel.h"
#include <fstream>

extern "C" {
        #include "dw_math.h"
        #include "dw_switchio.h"
        #include "dw_rand.h"
}

using namespace std; 

int CEquiEnergyModel::EE_Draw_RandomBlock(const CEESParameter &parameter, CStorageHead &storage)
{
        CSampleIDWeight x_new;
        int new_sample_code = NO_JUMP;
        double bounded_log_posterior_new;

        if (energy_level == parameter.number_energy_level-1 || dw_uniform_rnd() > parameter.pee)
        {
                if (metropolis->RandomBlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample))
                {
                        current_sample = x_new;
                        current_sample.id = (int)(time(NULL)-timer_when_started);
                        new_sample_code = METROPOLIS_JUMP;
                }
        }
	else
        {
		int bin_id = parameter.BinIndex(-current_sample.weight, energy_level+1);
                if (storage.DrawSample(bin_id, x_new)) // if a sample is successfully draw from bin
                {
			double log_ratio = parameter.LogRatio_Level(-x_new.weight, -current_sample.weight, energy_level);
                        log_ratio += parameter.LogRatio_Level(-current_sample.weight, -x_new.weight, energy_level+1);
			if (log(dw_uniform_rnd()) <= log_ratio)
                        {
				current_sample = x_new;
                                current_sample.id = (int)(time(NULL)-timer_when_started);
                                new_sample_code = EQUI_ENERGY_JUMP;
                        }
                        else
                        {
                                if (metropolis->RandomBlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample))
                                {
                                        current_sample = x_new;
                                        current_sample.id = (int)(time(NULL)-timer_when_started);
                                        new_sample_code = METROPOLIS_JUMP;
                                }
                        }
                }
		else
                {
                        if (metropolis->RandomBlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample) )
                        {
                                current_sample = x_new;
                                current_sample.id = (int)(time(NULL)-timer_when_started);
                                new_sample_code = METROPOLIS_JUMP;
                        }
                }
        }

        return new_sample_code;

}	

double CEquiEnergyModel::Simulation_Within_RandomBlock(const CEESParameter &parameter, CStorageHead &storage, bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new;
        unsigned int nJump =0;
        double max_posterior = current_sample.weight, bounded_log_posterior_new;
        bool if_write_file = false;
        ofstream output_file;
        if (!sample_file_name.empty() )
        {
                output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
                if (output_file)
                        if_write_file = true;
        }

        for (unsigned int i=0; i<parameter.simulation_length; i++)
        {
                for (unsigned int j=0; j<parameter.deposit_frequency; j++)
		if (metropolis->RandomBlockRandomWalkMetropolis(bounded_log_posterior_new, x_new, current_sample) )
        	{
			current_sample = x_new;
               		current_sample.id = (int)(time(NULL)-timer_when_started);
               		max_posterior = current_sample.weight > max_posterior ? current_sample.weight : max_posterior;
               		nJump ++;
        	}
		if (if_storage)
       			SaveSampleToStorage(current_sample, parameter, storage);
       		if (if_write_file)
      			 write(output_file, &current_sample);
       }
       if (if_write_file)
       			output_file.close();
	cout << "MH Jump " << nJump << " out of " << parameter.simulation_length << " in simulation.\n"; 
	return max_posterior;
}

double CEquiEnergyModel::Simulation_Cross_RandomBlock(const CEESParameter &parameter, CStorageHead &storage, bool if_storage, const string &sample_file_name)
{
	CSampleIDWeight x_new;

        unsigned int nEEJump=0, nMHJump=0;
        bool if_write_file = false;
        ofstream output_file;
        if (!sample_file_name.empty() )
        {
                output_file.open(sample_file_name.c_str(), ios::binary | ios::out);
                if (output_file)
                        if_write_file = true;
        }

        double max_posterior = current_sample.weight;
        for (unsigned int i=0; i<parameter.simulation_length; i++)
        {
                for (unsigned int j=0; j<parameter.deposit_frequency; j++)
                {
                        int jump_code = EE_Draw_RandomBlock(parameter, storage);
                        if (jump_code == EQUI_ENERGY_JUMP)
                                nEEJump++;
                        else if (jump_code == METROPOLIS_JUMP)
                                nMHJump++;
                        if (jump_code == EQUI_ENERGY_JUMP || jump_code == METROPOLIS_JUMP)
                                max_posterior = max_posterior > current_sample.weight ? max_posterior : current_sample.weight;
                }

                if (if_storage)
                        SaveSampleToStorage(current_sample, parameter, storage);
                if (if_write_file)
                        write(output_file, &current_sample);
        }
	if (if_write_file)
                output_file.close();

        cout << "EE Jump " << nEEJump << " out of " << parameter.simulation_length << " in simulation.\n";
        cout << "MH Jump " << nMHJump << " out of " << parameter.simulation_length << " in simulation.\n";
        return max_posterior;
}
