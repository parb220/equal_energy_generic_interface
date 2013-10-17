#include "CEquiEnergy_TState.h"
using namespace std; 

bool CEquiEnergy_TState::SaveTargetModelOriginalSetting()
{
        if (target_model == NULL)
                return false;

        size_t n = NumberFreeParametersTheta(target_model)+NumberFreeParametersQ(target_model);
        original_sample.data.Resize(n);
        double *x = new double[n];
        ConvertThetaToFreeParameters(target_model,x);
        ConvertQToFreeParameters(target_model,x+NumberFreeParametersTheta(target_model) );

        for (int i=0; i<original_sample.data.dim; i++)
                original_sample.data[i] = x[i];
        original_sample.DataChanged();
        original_sample.id = (int)(time(NULL)-timer_when_started);
        log_posterior_function(original_sample);
        return true;
}

bool CEquiEnergy_TState::RecoverTargetModelOriginalSetting()
{
        if (target_model == NULL)
                return false;
        ConvertFreeParametersToTheta(target_model, original_sample.data.vector);
        ConvertFreeParametersToQ(target_model, original_sample.data.vector+NumberFreeParametersTheta(target_model) );
        return true;
}

bool CEquiEnergy_TState::InitializeFromTarget()
{
        if (target_model == NULL)
                return false;
	// size_t n = NumberFreeParametersTheta(target_model)+NumberFreeParametersQ(target_model);
	// current_sample.data.Resize(n);
	// double *x = new double[n];
	// ConvertThetaToFreeParameters(target_model,x);
	// ConvertQToFreeParameters(target_model,x+NumberFreeParametersTheta(target_model) );
	
	// for (int i=0; i<current_sample.data.dim; i++)
	// current_sample.data[i] = x[i]; 
	// current_sample.DataChanged(); 
	current_sample = original_sample;
        current_sample.id = (int)(time(NULL)-timer_when_started);
	// log_posterior_function(current_sample);
	// delete []x;
	return true; 
}

double CEquiEnergy_TState::log_posterior_function(CSampleIDWeight &x)
{
        if (!x.calculated)
        {
		// double *old_x = new double[x.data.dim];
		// Save the old parameters stored in target_model to old_x
		// ConvertThetaToFreeParameters(target_model, old_x);
		// ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) );
		// Post x to target_model
		ConvertFreeParametersToTheta(target_model, x.data.vector);
                ConvertFreeParametersToQ(target_model, x.data.vector+NumberFreeParametersTheta(target_model) );
                x.weight = LogPosterior_StatesIntegratedOut(target_model);
                x.calculated = true;
		// Post old_x back to target_model 
		// ConvertFreeParametersToTheta(target_model, old_x); 
		// ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) );
		// delete [] old_x;
	}
		double bounded_log_posterior;
        if (if_bounded)
                bounded_log_posterior = x.weight/t_bound;
        else
                bounded_log_posterior = x.weight;

        return bounded_log_posterior;
}

double CEquiEnergy_TState::log_likelihood_function(const CSampleIDWeight &x)
{
	// double *old_x = new double[x.data.dim]; 
	// Save the old parameters stored in target_model to old_x
	// ConvertThetaToFreeParameters(target_model, old_x);
	// ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) );
	
	// post x to target_model
	ConvertFreeParametersToTheta(target_model, x.data.vector);
        ConvertFreeParametersToQ(target_model, x.data.vector+NumberFreeParametersTheta(target_model) );
        double log_likelihood = LogLikelihood_StatesIntegratedOut(target_model);

	// post old_x back to target_model
	// ConvertFreeParametersToTheta(target_model, old_x);
	// ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) );
	// delete [] old_x;
	
	return log_likelihood;
}

double CEquiEnergy_TState::log_posterior_function(const double *x, size_t n)
{
	//double *old_x = new double[n];
	//ConvertThetaToFreeParameters(target_model, old_x);
	//ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) );
	
	ConvertFreeParametersToTheta(target_model, (double*)x);
        ConvertFreeParametersToQ(target_model, (double*)x+NumberFreeParametersTheta(target_model) );
        double log_posterior = LogPosterior_StatesIntegratedOut(target_model);

	//ConvertFreeParametersToTheta(target_model, old_x);
	//ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) );
	//delete [] old_x;
	
	return if_bounded ? log_posterior/t_bound : log_posterior;
}

double CEquiEnergy_TState::log_likelihood_function(const double *x, size_t n)
{
	// double *old_x = new double[n];
	// ConvertThetaToFreeParameters(target_model, old_x);
	// ConvertQToFreeParameters(target_model, old_x+NumberFreeParametersTheta(target_model) );
	
	ConvertFreeParametersToTheta(target_model, (double *)x);
        ConvertFreeParametersToQ(target_model, (double *)x+NumberFreeParametersTheta(target_model) );
        double log_likelihood = LogLikelihood_StatesIntegratedOut(target_model);

	// post old_x back to target_model
	// ConvertFreeParametersToTheta(target_model, old_x);
	// ConvertFreeParametersToQ(target_model, old_x+NumberFreeParametersTheta(target_model) );
	// delete [] old_x;
	return log_likelihood;
}


CEquiEnergy_TState::CEquiEnergy_TState() : CEquiEnergyModel(), original_sample(CSampleIDWeight()), target_model(NULL)
{}

CEquiEnergy_TState::CEquiEnergy_TState(bool _if_bounded, int eL, double _t, const CSampleIDWeight &_x, time_t _time, CMetropolis *_metropolis, CEESParameter *_parameter, CStorageHead *_storage, TStateModel *_model) : CEquiEnergyModel(_if_bounded, eL, _t, _x, _time, _metropolis, _parameter, _storage), original_sample(CSampleIDWeight()), target_model(_model)
{
}

