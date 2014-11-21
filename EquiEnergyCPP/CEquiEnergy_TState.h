#ifndef _EQUI_ENERGY_TSTATE_
#define _EQUI_ENERGY_TSTATE_

#include "CEquiEnergyModel.h"
extern "C"
{
	#include "dw_switch.h"
}

class CEquiEnergy_TState : public CEquiEnergyModel
{
private:
	CEquiEnergy_TState(const CEquiEnergy_TState &); 
	const CEquiEnergy_TState & operator=(const CEquiEnergy_TState &); 
protected:
	CSampleIDWeight original_sample;
public: 
	TStateModel *target_model;
	bool SaveTargetModelOriginalSetting(); 
	bool RecoverTargetModelOriginalSetting();
	bool InitializeFromTarget(); 
	virtual double log_posterior_function(CSampleIDWeight &x); 
	// x cannot be constant because x.weight will be set as the real log_posterior calculated
	// from target_model, where the returning value is the bounded log_posterior
	// that is, return value = if_bounded ? x.weight/t_bound : x.weight; 
	virtual double log_likelihood_function(const CSampleIDWeight &x); 
	virtual double log_prior_function(const CSampleIDWeight &x); 
	// returning value is the real log_likelihood calculated from target_model
	virtual double log_posterior_function(const double *x, int n);
        virtual double log_likelihood_function(const double *x, int n);
        virtual double log_prior_function(const double *x, int n);

	using CEquiEnergyModel::HillClimb_NPSOL;
        // using CEquiEnergyModel:: HillClimb_CSMINWEL;
	virtual bool DrawParametersFromPrior(double *x) const; 

	// Construction 
	CEquiEnergy_TState(); 
	CEquiEnergy_TState(bool _if_bounded, int eL, double _t, const CSampleIDWeight &_x, time_t _time, CMetropolis *_metropolis, CEESParameter *_parameter, CStorageHead *_storage, TStateModel *_model); 
	~CEquiEnergy_TState() {}
};

#endif
