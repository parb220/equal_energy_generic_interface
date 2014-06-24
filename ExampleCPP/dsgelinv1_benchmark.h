#ifndef __DSGELINV1_SPECIFICMOD_H__
#define __DSGELINV1_SPECIFICMOD_H__

#include "dsgelinv1_comfuns.h"

//=== $$##:  All these parameters are model-specific.
typedef struct TSfinalpars_tag
{
   int n_constpars;  // number of constant (time-invariant) free parameters.
   //--- Deep constant parameters
   double gbeta;         //Pre-fixed value when indx_hard_restrict = 2 in datainp_common.prn.
   double glambdaq;      //Pre-fixed value when indx_hard_restrict = 207 in datainp_common.prn.
   double glambdastar;  //Pre-fixed value when indx_hard_restrict = 1,4 in datainp_common.prn.
   double ggammah;
   double ggammae;
   double geta;         //Pre-fixed value when indx_hard_restrict = 1,3,4 in datainp_common.prn.
   double gOmega;

   //--- Steady state.
   double Iql;
   double Iqk;
   double gthetaK;
   //---
   double R;         //Free parameter for mod5.
   double by;
   double ky;
   double ik;
   double qlLeOY;
   double qlLhOY;
   double nbar;

   //------ Derived parameter values. Will be imported from Cinput_DerivedPars.dat
   double glambdaK;
   double glambdaa;
   double gtheta  ;
   double gdelta  ;
   double iy      ;
   double galpha  ;     //Pre-fixed value when indx_hard_restrict = 3,4,5 in datainp_common.prn.
   double gphi    ;
   double CeOY    ;
   double ChOY    ;
   double gvarphi ;
   double glambdaz;
   //--- Inputs from Mathematica
   double gmuboe ;
   double qlLeOB ;
   double LhOLe  ;
   double Le     ;
   double Lh     ;
   double gOmegah;
   double gOmegae;
   double ER     ;       //Excess return.
   //--- A few more for computing the LH in Refresh_kalfilms_allcases() and in Cinput_MeasureEqns.exps.
   double InvdenTerm;
   double numTerm   ;
   double COY       ;
   double CeOC      ;
   double ChOC      ;

   //--- Exogenous persistence coefficients - time varying.
   int n_tvcoefpars; //number of coefficient parameters that are time-varying.
   TSdvector *grhoglambdaa_dv;
   TSdvector *grhoglambdaz_dv;
   TSdvector *grhognuz_dv;
   TSdvector *grhoglambdaq_dv;
   TSdvector *grhognuq_dv;
   TSdvector *grhogvarphi_dv;
   TSdvector *grhogpsi_dv;
   TSdvector *grhogtheta_dv;
   TSdvector *grhogdelta_dv;
   //--- Transition matrix -- not used for this project.
   TSdmatrix *Qc_dm;    //n_st_coef[0]-by-n_st_coef[0]: will be copied from smodel_ps->sv->QA[0];  //Transition matrix for coefficients changing.



   //--- Shock variances -- time varying.
   int n_tvvarpars; //number of variance parameters that are time varying (= nshocks).
   TSdvector *gsigmaglambdaa_dv;
   TSdvector *gsigmaglambdaz_dv;
   TSdvector *gsigmagnuz_dv;
   TSdvector *gsigmaglambdaq_dv;
   TSdvector *gsigmagnuq_dv;
   TSdvector *gsigmagvarphi_dv;
   TSdvector *gsigmagpsi_dv;
   TSdvector *gsigmagtheta_dv;
   TSdvector *gsigmagdelta_dv;
} TSfinalpars;


typedef struct TSpriorhyperpars_tag
{
   //======= Hyperparameters for deep parameters. =======
   double a_x_gbeta      ;
   double a_x_glambdaq   ;
   double a_x_glambdastar;
   double a_ggammah    ;
   double a_ggammae    ;
   double a_geta       ;
   double a_gOmega     ;
   double a_x_R        ;  //If Model 5 is used.
   //+
   double b_x_gbeta      ;
   double b_x_glambdaq   ;
   double b_x_glambdastar;
   double b_ggammah    ;
   double b_ggammae    ;
   double b_geta       ;
   double b_gOmega     ;
   double b_x_R        ;  //If Model 5 is used.


   //=== Exogenous processes
   // double a_grhoglambdaa;
   // double a_grhoglambdaz;
   // double a_grhognuz;
   // double a_grhoglambdaq;
   // double a_grhognuq;
   // double a_grhogvarphi;
   // double a_grhogpsi;
   // double a_grhogtheta;
   // //+
   // double b_grhoglambdaa;
   // double b_grhoglambdaz;
   // double b_grhognuz;
   // double b_grhoglambdaq;
   // double b_grhognuq;
   // double b_grhogvarphi;
   // double b_grhogpsi;
   // double b_grhogtheta;
   //
   double a_rho_all;
   double b_rho_all;


   //--- Shock deviatons: inverse-gamma prior for standard deviations of all shocks:
   double a_sigma_all;    // for all
   double b_sigma_all;    //Gamma prior

} TSpriorhyperpars;

//+
struct TSfinalpars_tag *CreateTSfinalpars(FILE *fptr_common, struct TSdsgemodel_tag *dsgemodel_ps);    //Called by dsgelinv1_estmcmc.c.
struct TSfinalpars_tag *DestroyTSfinalpars(struct TSfinalpars_tag *finalpars_ps);                      //Called by dsgelinv1_comfuns.c.
struct TSpriorhyperpars_tag *CreateTSpriorhyperpars(void);                                             //Called by dsgelinv1_comfuns.c.
////struct TSpriorhyperpars_tag *DestroyTSpriorhyperpars(struct TSpriorhyperpars_tag *priorhyperpars_ps);  //Called by dsgelinv1_comfuns.c.
double logpriordensity(struct TStateModel_tag *smodel_ps);                                             //Called by dsgelinv1_comfuns.c.
void Draw_xphi_dvFromPrior_Specific(struct TStateModel_tag *smodel_ps);                                //Called by dsgelinv1_estmcmc.c and dsgelinv1_comfuns.c.


//------------------------------------------------------------
// Functions required by dw_temp_output.c in DW's switching subdirectory
//   for DW's kalman-filtering functions.
//------------------------------------------------------------
void Convertphifreepars2finalpars(struct TSdsgemodel_tag *dsgemodel_ps);  //Converting free parameters (bar transition matrix elements) to the final (DSGE) meaningful parameters.
                                                                          //Called by dsgelinv1_comfuns.c.
int ViolateAPrioriRestrictions(struct TSdsgemodel_tag *dsgemodel_ps);     //Called by dsgelinv1_comfuns.c.



#endif
