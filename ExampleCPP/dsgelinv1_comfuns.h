//-------------- Common Functions ---------------------
//-- Looking for $$## to search for hard-coded items.
//--- In CreateTSpriorhyperpars(), Prior 2 is new and Prior 1 corresponds to old prior (with less persistence).

#ifndef __DSGELINV1_COMFUNS_H__
#define __DSGELINV1_COMFUNS_H__

#include "tzmatlab.h"
#include "fn_filesetup.h"
#include "mathlib.h"
#include "gensys.h"
#include "kalman.h"
#include "cstz.h"
#include "rand.h"
#include "optpackage.h"
#include "switch.h"
#include "switchio.h"
#include "dw_array.h"
//--- Added by TZ
#include "dw_metropolis_simulation.h"
#include "dw_state_space_impulse_response.h"
#include "dw_mdd_switch.h"
#include "dw_MSStateSpace.h"
#include "dw_state_space_forecast.h"  //Used only for dw_state_space_forecast_command_line() by dsgelinv1lwz_estmcmc.c --  9 April 2010.
#include "dw_parse_cmd.h"   //All commandline functions such as dw_ParseString_String(n_args_cl, args_cl, "random",(char*)NULL)
//--- The order in the following three lines matters!  That is, these "struct"s are put in front of #include "dsgelinv1_cexps.h"
//---   because this .h file needs to know these structs.
struct TSdsgemodel_tag;
struct TSDWRegimeControlVariables_tag;
#include "dsgelinv1_cexps.h"



//#define DEBUG_ON    //No longer used, 12 April 2010.  Replaced by FPTR_DEBUG in the main.c file.  //If define, in debug mode; else (undef), in a normal running mode.
#define TAGLEN 512


//************************************************************
// Grand structure for the DSGE model.  p7d: parameters and data.
//************************************************************
typedef struct TSdsgemodel_tag
{
   //=== From the data file dataraw*.prn.
   int ny;   //Number of observables in CreateTSkalfilmsinputs_1stapp2() in kalman.c.
   int nz;   //Number of state variables in CreateTSkalfilmsinputs_1stapp2() in kalman.c.  May > ngensys because of non-stationary variables.
   int nzbase;  //Number of base state variables in Kalman filtering (i.e., no attached lagged variables or without time-varying inflation target (pi*(s_t)) and regime (ebar(s_t)) variables).
   int nadded; //nadded = nz-nzbase;  nzbase <= ngensys but nzbase = ngensys if no regime-swtiching constant terms in the gnesys system.
   //int ntrendvars; //$$## number of trend (non-stationary) variables used in measurement equations.  MANUALLy changed with model.
   //int naddvars_kalman; //$$## number of additional variables in state equations used for Kalman filtering.  See f_t on p.48.
   //int nbasestatevars; //$$## number of base state variables in Kalman filtering (i.e., without time-varying inflation target (pi*(s_t)) and regime (ebar(s_t)) variables).
                    // naddvars_kalman + nbasestatevars = nz (total number in the Kalman filtering).
   int nu;   //number of measurement errors in Kalman filtering in CreateTSkalfilmsinputs_1stapp2() in kalman.c.
   int ne;   //number of fundamental shocks in Kalman filtering in CreateTSkalfilmsinputs_1stapp2() in kalman.c.
             //  May > nshocks in gensys when some AR measurement errors are grouped as state equations.

   int nDataVars;  //number of variables in the raw data.
   int nData;  //nnans+nlags+fss.
   int nnans;  //number of NaNs or data points that are NOT used.
   int nlags;  //total number of lags where the kalman filtering go through nlags and the sample starts at nlags+1.
   int fss;    //effective sample size = nData-nnans-nlags.
   int nSample;  //nlags+fss.
   //--- Reads from dataraw_*_logData_all.prn.
   TSdmatrix *Dataraw_dm;  //nData-by-nDataVars.
   TSdmatrix *Datarawtran_dm;  //nDataVars-by-nData. All the data we have.  
   TSivector *BeginEndData_iv;  //4-by-1: yearBegin monthBegin yearEnd monthEnd  All the data we have.
   //--- Reads from datainpu_TAG.prn.
   TSdmatrix *SampleData_dm;  //(fss+nlags)-by-n_y where n_y is number of observables.
   TSdmatrix *SampleDatatran_dm;  //n_y-by-(fss+nlags) where n_y is number of observables.
   TSivector *BeginEndSample_iv;  //4-by-1: yearBegin monthBegin yearEnd monthEnd.  Specific sample chosen for estimation.  The length is fss+nlags.
   TSivector *SelectData_iv;   //Selecting the appropriate data according to dataraw_*_logData_all.prn.  Must match n_y.
   TSdvector *dates_dv;  //(fss+nlags)-by-1: dates for the whole sample (including lags).



   //---------------------- Specific model ---------------------------------
   TScvector *modeltag_cv;  //allocated to 512 for a possible long name or tag.
   int indx_hard_restrict; //0: default model (persistence for all exogenous shocks)
                           //1: restricting glambdastar and geta but everything is the same as model 0 (default model).
                           //2:
   int flag_ComputingHessian;  //1: computing Hessian; 0: no computing.
   int flag_ComputingHist;     //1: Computing historical decompositions (and forecasts); 0: No computing when doing estimation.
   int indx_tvmodel; //0: no time-varying coefficient parameters present (only variance parameters may be time-varying).
                     //1: there are some time-varying coefficient parameters present.
   int indx_const;   //0: some time-varying parameters (either coefficients or shock variances).
                     //1: no regime-switching at all (not even regime-switching shock variances).
   int refreshed;   //IMPORTANT: this value will be changed to 0 WHENEVER tz_thetaChanged() or
                    //   tz_TransitionMatrixChanged() is called.  These functions are most likely to be called
                    //   by DW's internally consistent code.
                    //1: everything that needs be computed given new parameters has already been computed
                    //   (no need to recompute. For example, computing gensys once while summing over
                    //    the logLH_timet over t, so that gensys won't be computed every time).
                    //0: must compute everything that needs be computed given new parameters.
   int ValidParameters;  //1: parameters (bar transition matrix elements) are valid (satisfying a priori restrictions),
                         //0: not valid.
                         //This has a meaning ONLY when refreshed=1.
   int ValidGensysSolution;  //1: gensys has a (MSV or unique) solution,
                             //0: no such a solution,
                             //This has a meaning ONLY when refreshed=1.
   int ValidInitial_z10_P10; //1: InitializeKalman_z10_P10() gives valid initial values z_{1|0} and P_{1|0}.
                             //0: no such valid initials.
                             //This has a meaning ONLY when refreshed=1.


   //--- Minimization problem.
   int indxStartValuesForMin;  //1: initial start from datainp_setup.prn; 0: continuing from the last estimated results if they exist.
   int indxEstFinder;   //1: finds the posterior or ML estimate; 0: no such a finding and imports values of the parameters from a file.
   int nDrawsInitsp;  //Number of initial starting points draws randomly from the prior to find the posterior peak (with the option /c 2).
   ////int indxSimulateMCMC;  //1: simulating MCMC using the DW code (can be time-consuming); 0: no such simulations.
   double LHvalatpostpeak;  //LH value at the posterior peak.
   double peaklogpost;
   double logprior;
   //
   int randomseed;
   time_t prog_begtime;
   time_t prog_endtime;
   double program_hours;      //Hours spent on the program.

   //--- Drawing from the prior.
   int nDrawsFromPrior; //Number of draws from the prior to determine the scale.
   int nViolates;  //Number of violations of the restrictions from the draws from the prior.
   double scale4logpriordensity;  //Scale of log prior density for possible hard restrictions on the parameter space.

   //=== DW regime control variables.
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps;

   //=== Final parameters used by the model.
   struct TSfinalpars_tag *finalpars_ps;  //free parameters (time-invariant and time-varying).
   //--- Read-in initial free parameters used for finding maximization, allocataed memory only if indxStartValuesForMin==1.
   TSdvector *xphi_dv;   //free parameters (NOT including transition matrix elements) for MCMC; its dimension (depending on indx_tvmodel) will always be read.  For the vector itself, read in only if indxStartValuesForMin==1.
   //--- Covariance matrix for standard errors at the log posterior peak.
   int nfreepars_tot; //A total number of free parameters (length(xphi_dv)+number of free Q).
   TSdmatrix *Omega_dm;   //nfreepars_tot-by-nfreepars_tot or minpack_ps->x_dv->n-by-minpack_ps->x_dv->n.
   //--- Hessian matrix at the log posterior peak -- its inverse = Omega (covariance matrix).
   TSdmatrix *Hessian_dm;   //nfreepars_tot-by-nfreepars_tot or minpack_ps->x_dv->n-by-minpack_ps->x_dv->n.

   //--- Gensys structure.
   int ngensys;  //Total number of variables in the gensys form.  ngensys may NOT equal nz in tz_kalfiltv() in kalman.c because lagged variables may be added to the system after the gensys is solved.
                   //  ngensys may not equal nbasestatevars because regime variables may add to the system for gensys and then later move to residuals after the gensys is solved.
   int nshocks;    //Number of all i.i.d. shocks (u_t) for the gensys form: number of fundamental shocks + number of coefficient regimes.
   int nshocks_kalman;  //Number of fundamental shocks for the final state space form for Kalman filter.
   int nexperrs;  //Number of expectational errors (eta_t) in gensys.
   double div;   //Division for separating large vs small roots.  Used only for gensys_sims.
   struct TSgensys_tag *gensys_ps;


   //=== IMSL linear constraints.
   int Use_IMSL_opt;  //1: use the IMSL optimization package; 0: not use the package so the program will NOT read all the following settings.
   int imsl_neqs; //number of equality constraints, excluding simple bound constraints.  IMSL dictates that equality constraints come always BEFORE inequality constrains.
   int imsl_ncons; //a total number of constrains, including equality and inequality constraints, but excluding simple bound constraints.
   TSivector *imsl_lh_locs_iv; //base-0 locations out of the npars*ncons vector for the nonzero coefficients in left-hand sides of the linear constraints.
   TSdvector *imsl_lh_vals_dv; //values of the nonzero coefficients corresponding to the locations specified by imsl_lh_locs_dv  in left-hand sides of the linear constraints.
   TSdvector *imsl_rh_vals_dv; //values in right-hand sides of the linear constraints. The number of these values matches imsl_ncons.
   //+
   int hard7_imsl_neqs; //number of equality constraints, excluding simple bound constraints.  IMSL dictates that equality constraints come always BEFORE inequality constrains.
   int hard7_imsl_ncons; //a total number of constrains, including equality and inequality constraints, but excluding simple bound constraints.
   TSivector *hard7_imsl_lh_locs_iv; //base-0 locations out of the npars*ncons vector for the nonzero coefficients in left-hand sides of the linear constraints.
   TSdvector *hard7_imsl_lh_vals_dv; //values of the nonzero coefficients corresponding to the locations specified by hard7_imsl_lh_locs_dv  in left-hand sides of the linear constraints.
   TSdvector *hard7_imsl_rh_vals_dv; //values in right-hand sides of the linear constraints. The number of these values matches hard7_imsl_ncons.
   //--- Simple bounds for the default model (persistence for all exogenous shocks)
   TSivector *hard0_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard0_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard0_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for model 1
   TSivector *hard1_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard1_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard1_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for model 2
   TSivector *hard2_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard2_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard2_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for model 3
   TSivector *hard3_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard3_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard3_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for modelb 103
   TSivector *hard103_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard103_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard103_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for modelb 203
   TSivector *hard203_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard203_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard203_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for model 4
   TSivector *hard4_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard4_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard4_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for model 5
   TSivector *hard5_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard5_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard5_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for modelb 105
   TSivector *hard105_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard105_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard105_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   //--- Simple bounds for modelb 205
   TSivector *hard205_simple_locs_iv; //base-0 locations out of the npars vector for simple-bound constraints.
   TSdvector *hard205_simple_lowvals_dv; //low values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.
   TSdvector *hard205_simple_highvals_dv; //high values of simple constraints, corresponding to the locations specified by hard#_simple_locs_iv.

  
   //=== NPSOL linear constraints.
   int Use_NPSOL_opt;  // 1: use the NPSOL optimization package; 0: not use the package so the program will NOT read all the following settings.
   int npsol_nlin;     // The number of linear constraints
   int npsol_nclin;    // The number of nonlinear constraints
   //--- New Constraint Structure
   struct TSModelConstraints_tag *model_constraints;



   //--- Markov-switching Kalman filter structure.
   struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps; //For getting initial conditions using the lagged data.

   //--- Structure for hyperparameters for the prior of model free parameters.
   struct TSpriorhyperpars_tag *priorhyperpars_ps;


   //=== Regime-switching output arguments.
   TSdvector *dw_ProbStates_dv;  //nData+1-by-1, calling DW's ProbabilitiesState() and then copying this vector to PostProbS_dm in appropriate places.
   TSdmatrix *PostProbS_dm;  //(fss+1)-by-smodel_ps->sv->nstates: marginal posterior probabilities of states starting
                             //  with [1] (base-1) and with [0] as an initial state.
   TSdcell *BaseProbS_dc;  //(fss+1)-by-nbasestates-by-nstatevariables: marginal posterior probabilities of states starting
                             //  with [1] (base-1) and with [0] as an initial state.

   //--- Not really right for historical decompositon or computing the structural shocks.
   TSdmatrix *strshks_dm; //ny-by-fss, structural shocks assuming no measurement errors.     ?????? Questionable historical decomps based on Junior's notes.
   // TSdcell *Impact_t_exact_dc; //nz-by-ny-by-nst.  Impact with exactly ny fundamental shocks.  ?????? Questionable historical decomps based on Junior's notes.

} TSdsgemodel;


typedef struct TSDWRegimeControlVariables_tag
{
   //--- Memory of all the following arrays will be allocated by DW and destroyed by the DW function dw_FreeArray(),
   //---   partly because DW allocated array with some information stored in the -2, -1 elements, so
   //---   the array does not begin with base 0.   Thus, tzDestroy(n_coef) or free(n_coef)will NOT work.
   //--- Thus, it's essential to destroy them using dw_FreeArray().

   int n_constpars;  // number of constant (time-invariant) free parameters.

   int *n_st_coef;   //n_st_coef[i-1] is the number of regimes associated with the i_th time-varying coefficient parameter.
   int *cum_total_coef; //cum_total_coef[i-1] is the offset of the i_th time-varying coefficient parameter
                        //  in the whole vector of free parameters.
   int **indx_st_coef; //indx_st_coef[st][i-1] is the regime associated with the i_th time-varying
                       //  variance parameter with the grand regime st.

   int *n_st_var; //n_st_var[i-1] is the number of regimes associated with the i_th time-varying variance parameter.
   int *cum_total_var; //cum_total_var[i-1] is the offset of the i_th time-varying variance parameter
                       //  in the whole vector of free parameters.
   int **indx_st_var; //indx_var[st][i-1] is the regime associated with the i_th time-varying
                      //  variance parameter with the grand regime st.
                     //Use 1: g_sigma_w_dv->v[indx_st_var[st][2]] where 2 means that g_sigma_w_dv is the 3rd
                     //         time-varying variance parameter to be considered.
                     //Use 2: xphi_dv->v[cum_total_var[2]+indx_st_var[st][2]] where 2 means the 3rd time-varying
                     //         variance parameters;  this value is the same as g_sigma_w_dv->v[idx_var[st][2]].
} TSDWRegimeControlVariables;





//--- Declaration of protypes.
struct TSdsgemodel_tag *CreateTSdsgemodel(FILE *fptr_rawdata, FILE *fptr_common, FILE *fptr_tvpars, FILE *fptr_markov, int indxStartValuesForMin);
struct TSdsgemodel_tag *DestroyTSdsgemodel(struct TSdsgemodel_tag *dsgemodel_ps);
void InitializeCreateTSdsgemodelFromStateModel(struct TStateModel_tag *smodel_ps);
struct TSDWRegimeControlVariables_tag *CreateTSDWRegimeControlVariables(FILE *fptr_markov, struct TStateModel_tag *smodel_ps);
void FindFatalErrorsOnDimensions(struct TSdsgemodel_tag *dsgemodel_ps);
//
int RefreshEverything(struct TStateModel_tag *smodel_ps);
double tz_GetScaleForLogpriordensity(struct TStateModel_tag *smodel_ps);
double logTimetCondLH(int s_t, int inpt, struct TStateModel_tag *smodel_ps);
double logOverallPosteriorKernal_const(struct TStateModel_tag *smodel_ps, TSdvector *xnew_dv);


//--- Free Parameter Structure
struct TSFreeParameter_tag {
  char *name;
  char *variable_name;
  char *prior_distribution;
  TSdvector *prior_parameters;
  float  lower_bound;
  float  upper_bound;
} TSFreeParameter;

struct TSFreeParameter_tag *CreateTSFreeParameter(char *name);

//--- New Constraint Types
struct TSModelConstraints_tag {
  int n;      // The number of parameters
  int num_lin;   // Number of linear constraints
  int num_nlin;  // The number of nonlinear constraints
  
  // General bounds on parameters
  TSdvector *lb;
  TSdvector *ub;
  
  // Linear Constraint Matrix
  TSdmatrix *A;
  // linear constraints lower_bounds
  TSdvector *lin_lb;
  TSdvector *lin_ub;
  
  // Non Linear constraints bounds
  TSdvector *nlin_lb;
  TSdvector *nlin_ub;
  
  void (*non_lin_const_func)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate);
  
} TSModelConstraints;
typedef void TFNonLinConstraints (int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate);
struct TSModelConstraints_tag *CreateTSModelConstraints(int n, int nlin, int nclin, TFNonLinConstraints *non_lin_con_func);



//------------------------------------------------------------
// Functions required by CreateParameters() in switch.c.
//------------------------------------------------------------
int NumberOfFreeModelSpecificParameters(struct TStateModel_tag *smodel_ps);
void ConvertFreeParameters2ModelSpecificParameters(struct TStateModel_tag *smodel_ps, PRECISION *xopt_pd);
void CopyMyFreeParameters2OptimizationParameters(struct TStateModel_tag *smodel_ps, PRECISION *xopt_pd);
void tz_thetaChanged(struct TStateModel_tag *smodel_ps);
void tz_TransitionMatrixChanged(struct TStateModel_tag *smodel_ps);







/***************************************************************************************/
//-------------------               IMPORTANT NOTES                            --------------------
//------------------  Steps for coding up regime-switching models: an example. --------------------

//1.   Creatae datainp_trendinf_2c2v.prn and make sure all numbers are correct for TSdsgemodel_tag.
//2.   Create TSfinalpars_tag, TSpriorhyperpars_tag, and TSdsgemodel_tag
//              (which creates CreateTSgensys, kalfilmsinputs_1stapp_ps, etc.).
//3    Change FindFatalErrorsOnDimensions();
//4.   Create RunningGensys_allcases().
//5.   Change Convertphifreepars2finalpars().
//6.   Change ViolateAPrioriRestrictions().
//7.   Create Refresh_kalfilms_allcases(), paying attention to c(s_{t-1}, c_t) or at_dm.
//8.   Change logpriordensity().
//9.   Examine to make sure the library tz_logTimetCondLH_kalfilms_1st_approx() in kalman.c is still good. ?????
//10.  Change RefreshEverything().  This is a very important function.  Check carefully.
//12.  Change logTimetCondLH().
//13.  Modify ftd_FprintOutput() and check indx_tvmodel.

/****************************************************************************************/

#endif



