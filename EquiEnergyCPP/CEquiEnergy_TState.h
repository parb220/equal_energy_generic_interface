#ifndef _EQUI_ENERGY_TSTATE_
#define _EQUI_ENERGY_TSTATE_

#include "CEquiEnergyModel.h"
extern "C"
{
	#include "dw_switch.h"
}

class MinusLogPosterior_NPSOL; 
class MinusLogPosterior_CSMINWEL; 
class CEquiEnergy_TState; 

class MinusLogPosterior_NPSOL
{
public:
        static CEquiEnergy_TState *model;
        static void *function(int *mode, int *n, double *x, double *f, double *g, int *nstate);
};

class MinusLogPosterior_CSMINWEL
{
public:
        static CEquiEnergy_TState *model;
        static double function(double *x, int n, double **args, int *dims);
};


class CEquiEnergy_TState : public CEquiEnergyModel
{
private:
	CEquiEnergy_TState(const CEquiEnergy_TState &); 
	const CEquiEnergy_TState & operator=(const CEquiEnergy_TState &); 
protected:
        double log_posterior_function(const double *x, size_t n);
        double log_likelihood_function(const double *x, size_t n);
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
	// returning value is the real log_likelihood calculated from target_model

	double HillClimb_NPSOL(size_t nSolution);
        double HillClimb_CSMINWEL(size_t nSolution) ;
	CEquiEnergy_TState(); 
	CEquiEnergy_TState(bool _if_bounded, int eL, double _t, const CSampleIDWeight &_x, time_t _time, CMetropolis *_metropolis, CEESParameter *_parameter, CStorageHead *_storage, TStateModel *_model); 
	~CEquiEnergy_TState() {}
friend class MinusLogPosterior_NPSOL;
friend class MinusLogPosterior_CSMINWEL;
};

#endif
