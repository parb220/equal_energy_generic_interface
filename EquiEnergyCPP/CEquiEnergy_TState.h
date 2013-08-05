#ifndef _EQUI_ENERGY_TSTATE_
#define _EQUI_ENERGY_TSTATE_

#include "CEquiEnergyModel.h"

class CEquiEnergy_TState : public CEquiEnergyModel
{
private:
	CEquiEnergy_TState(const CEquiEnergy_TState &); 
	const CEquiEnergy_TState & operator=(const CEquiEnergy_TState &); 
protected:
        virtual double log_posterior_function(const double *x, size_t n);
        virtual double log_likelihood_function(const double *x, size_t n);
	CSampleIDWeight original_sample; 
public: 
	TStateModel *target_model;
	virtual bool SaveTargetModelOriginalSetting(); 
	virtual bool RecoverTargetModelOriginalSetting();
	virtual bool InitializeFromTarget(); 
	virtual double log_posterior_function(CSampleIDWeight &x); 
	// x cannot be constant because x.weight will be set as the real log_posterior calculated
	// from target_model, where the returning value is the bounded log_posterior
	// that is, return value = -(-x.weight>h_bound ? -x.weight:h_bound)/t_bound;
	virtual double log_likelihood_function(const CSampleIDWeight &x); 
	// returning value is the real log_likelihood calculated from target_model

	CEquiEnergy_TState(); 
	CEquiEnergy_TState(bool _if_bounded, unsigned int eL, double _h, double _t, const CSampleIDWeight &_x, CMetropolis *_metropolis, time_t _time, TStateModel *_model); 
	~CEquiEnergy_TState() {}
};

#endif
