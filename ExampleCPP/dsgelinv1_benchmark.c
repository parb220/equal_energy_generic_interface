/*
 * Check $$## for hard-coded numbers or lines.
 * Check ??? for debugging or for unfinished jobs.
 * Check <<>> for updating DW's new switch code or questions for DW.
 *
 * Notes:
 *    indxStartValuesForMin: 1: initial start from datainp_setup.prn; 0: continuing from the last estimated results if they exist.
 *
 * ===== Steps for changing the model when indx_hard_restrict != 0 or when indx_hard_restrict is changed: ======
 *    1. Modify indx_hard_restrict, pre-fixed values, hard#_simple_locs_iv, etc in in datainp_common.prn
 *    2. Create datainp_markov_mod*_const.prn and datainp_mod*_const.prn (e.g, xphi_dv and Number of free parameters not time varying).
 *    3. Search for indx_hard_restrict in both *_comfuns.c and *_est.c.
 *    4. Modify the comments on parameters in TSfinalpars_tag in *_comfuns.h.
 *    5. Check if any variables to add (do nothing for taking off variables) in TSpriorhyperpars_tag in the *_comfuns.h file.
 *    6. Search for if any variables to add (do nothing for taking off variables) in CreateTSpriorhyperpars() in *_comfuns.c.
 *    7. Create gen_inipars_mod?_*.m
 *
 * ===== Ad hoc functions that must be changed for each new model or each different hard restrictions on parameters: ======
 *    Check in tzmatlab.h
 *          #define USE_DEBUG_FILE
 *          #define NEWVERSIONofDW_SWITCH
 *    Search for indx_hard_restrict in *_comfuns.c
 *    ----
 *    datainp_TAG.prn:
 *       //== xphi_dv ==//,
 *    datainp_markov_TAG.prn:
 *       //== Number of free parameters not time varying ==//
 *    datainp_common.prn:
 *        //== hard#_* ==//
 *    gen_inipars_TAG.m:
 *        When some values (e.g., R) in datainp_common.prn are changed, we need to change this value in gen_inipars_TAG.m
 *        accordingly.
 *    ----
 *    0.0. dsgelinv1_cexps.h.
 *    0.1. typedef struct TSfinalpars_tag in the .h file.
 *    0.2. CreateTSfinalpars() -- checking various //== ??? ==// read from datainp_common.prn.
 *    0.3. DestroyTSfinalpars()
 *    0.4. typedef struct TSpriorhyperpars_tag in the .h file.
 *    ----
 *    1. void Convertphifreepars2finalpars()
 *    2. void RunningGensys_allcases()
 *    3. void Refresh_kalfilms_allcases()   //Meaning: specifying measurement equations.
 *    4. CreateTSpriorhyperpars(void): checking the hyperparameter values to allow unit roots for technology, for example.
 *    5. logpriordensity()
 *    6. ViolateAPrioriRestrictions()
 *    7. ftd_FprintOutput() in the main C file.
 *    ----
 *    8. Looking for $$## to change hard-coded lines:
 *       FindFatalErrorsOnDimensions():
 *         if (finalpars_ps->n_constpars != 28)  //$$##
 *    9. Examine and modify
 *            Use_IMSL_opt, imsl_neqs, etc. and hard0_simple_locs_iv, etc
 *         in datainp_common.prn.
 *    10. Examine datainp_common.prn, datainp_markov_TAG.prn, and datainp_TAB.prn.
 *
 * Common functions for the linearized DSGE project.
 * Written by T. Zha, March 2009. Revisions: April 2010.
 * Copywright (c) by Zha, March 2009.
 *
*/

/**
//=== For debugging purpose.
if (1)
{
   double t_loglht;
   double jnk;

   jnk = logTimetCondLH_kalfilms_1stapp(s_t, inpt, kalfilmsinputs_1stapp_ps, smodel_ps);

   fprintf(FPTR_DEBUG, "\n-----------\n");
   fprintf(FPTR_DEBUG, "logLkelihood(s_t=%d, inpt=%d):\n", s_t, inpt);
   fprintf(FPTR_DEBUG, " %10.5f\n", jnk);

   fflush(FPTR_DEBUG);

   t_loglht = -(0.5*ny)*LOG2PI - 0.5*logdeterminant(Dtdata_dm) - 0.5*VectorDotVector(wny_dv, etdata_dv);
   fprintf(FPTR_DEBUG, " %10.5f\n", t_loglht);

   fprintf(FPTR_DEBUG, "%%st=%d, inpt=%d, and sti=%d\n", st, inpt, sti);

   fprintf(FPTR_DEBUG, "wP0_dv:\n");
   WriteVector(FPTR_DEBUG, wP0_dv, " %10.5f ");
   fprintf(FPTR_DEBUG, "Vt_dc->C[sti_v=%d]:\n", sti_v);
   WriteMatrix(FPTR_DEBUG, Vt_dc->C[sti_v], " %10.5f ");

   fflush(FPTR_DEBUG);
}
/**/

#include "dsgelinv1_comfuns.h"

struct TSfinalpars_tag *CreateTSfinalpars(FILE *fptr_common, struct TSdsgemodel_tag *dsgemodel_ps)
{
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   int *n_st_var = DWRegimeControlVariables_ps->n_st_var;
   //~~~ Creating the structure and initializing the NULL pointers.
   struct TSfinalpars_tag *finalpars_ps = tzMalloc(1,struct TSfinalpars_tag);

   if ( !fn_SetFilePosition(fptr_common, "//== Iql ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->Iql) != 1 )
      finalpars_ps->Iql = 1.0;
   if ( !fn_SetFilePosition(fptr_common, "//== Iqk ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->Iqk) != 1 )
      finalpars_ps->Iqk = 1.0;
   if ( !fn_SetFilePosition(fptr_common, "//== gthetaK ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->gthetaK) != 1 )
      finalpars_ps->gthetaK = 1.0;
   //---
   if ((dsgemodel_ps->indx_hard_restrict != 5) && (dsgemodel_ps->indx_hard_restrict != 105) && (dsgemodel_ps->indx_hard_restrict != 205)) //NOT Model 5.
      if ( !fn_SetFilePosition(fptr_common, "//== R ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->R) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to R");
   if (dsgemodel_ps->indx_hard_restrict <100)
      if ( !fn_SetFilePosition(fptr_common, "//== by ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->by) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to by");
   if ( !fn_SetFilePosition(fptr_common, "//== ky ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->ky) != 1 )
      fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to ky");
   if ( !fn_SetFilePosition(fptr_common, "//== ik ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->ik) != 1 )
      fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to ik");
   if ( !fn_SetFilePosition(fptr_common, "//== qlLeOY ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->qlLeOY) != 1 )
      fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to qlLeOY");
   if ( !fn_SetFilePosition(fptr_common, "//== qlLhOY ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->qlLhOY) != 1 )
      fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to qlLhOY");
   if ( !fn_SetFilePosition(fptr_common, "//== nbar ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->nbar) != 1 )
      fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to nbar");

   if (dsgemodel_ps->indx_hard_restrict==1)  //Model 1.
   {
      if ( !fn_SetFilePosition(fptr_common, "//== glambdastar ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->glambdastar) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to glambdastar");
      if ( !fn_SetFilePosition(fptr_common, "//== geta ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->geta) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to geta");
   }
   else if (dsgemodel_ps->indx_hard_restrict==2)  //Model 2.
   {
      if ( !fn_SetFilePosition(fptr_common, "//== gbeta ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->gbeta) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to gbeta");
   }
   else if (dsgemodel_ps->indx_hard_restrict==3 || dsgemodel_ps->indx_hard_restrict==5)  //Model 3 or 5.
   {
      if ( !fn_SetFilePosition(fptr_common, "//== galpha ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->galpha) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to galpha");
      if ( !fn_SetFilePosition(fptr_common, "//== geta ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->geta) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to geta");
   }
   else if (dsgemodel_ps->indx_hard_restrict==103 || dsgemodel_ps->indx_hard_restrict==105 ||
            dsgemodel_ps->indx_hard_restrict==203 || dsgemodel_ps->indx_hard_restrict==205)
   {
      if ( !fn_SetFilePosition(fptr_common, "//== galpha ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->galpha) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to galpha");
      if ( !fn_SetFilePosition(fptr_common, "//== geta ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->geta) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to geta");
      if ( !fn_SetFilePosition(fptr_common, "//== gtheta ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->gtheta) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to gtheta");
   }
   else if (dsgemodel_ps->indx_hard_restrict==207)
   {
      if ( !fn_SetFilePosition(fptr_common, "//== galpha ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->galpha) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to galpha");
      if ( !fn_SetFilePosition(fptr_common, "//== geta ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->geta) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to geta");
      if ( !fn_SetFilePosition(fptr_common, "//== gtheta ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->gtheta) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to gtheta");
      if ( !fn_SetFilePosition(fptr_common, "//== glambdaq ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->glambdaq) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to glambdaq");
      if ( !fn_SetFilePosition(fptr_common, "//== glambdastar ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->glambdastar) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to glambdastar");
   }
   else if (dsgemodel_ps->indx_hard_restrict==4)  //Model 4.
   {
      if ( !fn_SetFilePosition(fptr_common, "//== galpha ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->galpha) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to galpha");
      if ( !fn_SetFilePosition(fptr_common, "//== glambdastar ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->glambdastar) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to glambdastar");
      if ( !fn_SetFilePosition(fptr_common, "//== geta ==//") || fscanf(fptr_common, " %lf ", &finalpars_ps->geta) != 1 )
         fn_DisplayError(".../dsgelinv1_MODEL.c/CreateTSfinalpars(): need to assign a numer to geta");
   }



   //--- Number of time-invariant free parameters.
   finalpars_ps->n_constpars = DWRegimeControlVariables_ps->n_constpars;

   //=== Time-varying coefficients.
   if (dsgemodel_ps->indx_hard_restrict < 100)
      finalpars_ps->n_tvcoefpars = 8; //$$##: number of coefficient parameters that are time varying.
   else if (dsgemodel_ps->indx_hard_restrict < 200)
      finalpars_ps->n_tvcoefpars = 9; //$$##: number of coefficient parameters that are time varying.
   else if (dsgemodel_ps->indx_hard_restrict < 300)
      finalpars_ps->n_tvcoefpars = 8; //$$##: number of coefficient parameters that are time varying.
   //+
   if (dw_DimA(n_st_coef) != finalpars_ps->n_tvcoefpars)
      fn_DisplayError(".../*_comfuns.c/CreateTSfinalpars(): number of tv coefficients do not match");
   //+
   finalpars_ps->grhoglambdaa_dv = CreateConstantVector_lf(n_st_coef[0], 0.0);
   finalpars_ps->grhoglambdaz_dv = CreateConstantVector_lf(n_st_coef[1], 0.0);
   finalpars_ps->grhognuz_dv =     CreateConstantVector_lf(n_st_coef[2], 0.0);
   finalpars_ps->grhoglambdaq_dv = CreateConstantVector_lf(n_st_coef[3], 0.0);
   finalpars_ps->grhognuq_dv =     CreateConstantVector_lf(n_st_coef[4], 0.0);
   finalpars_ps->grhogvarphi_dv =  CreateConstantVector_lf(n_st_coef[5], 0.0);
   finalpars_ps->grhogpsi_dv =     CreateConstantVector_lf(n_st_coef[6], 0.0);
   finalpars_ps->grhogtheta_dv =   CreateConstantVector_lf(n_st_coef[7], 0.0);
   if ((dsgemodel_ps->indx_hard_restrict >= 100) && (dsgemodel_ps->indx_hard_restrict < 200))
      finalpars_ps->grhogdelta_dv =   CreateConstantVector_lf(n_st_coef[8], 0.0);
   else
      finalpars_ps->grhogdelta_dv = NULL;
                            //$$## hard-coded for [7] or [8].
   finalpars_ps->Qc_dm = NULL; //CreateMatrix_lf(n_st_coef[0], n_st_coef[0]);    //Not used.  Reserved for future use.
                            //Transition matrix


   //--- Shock standard deviations.
   if (dsgemodel_ps->indx_hard_restrict < 100)
      finalpars_ps->n_tvvarpars = 8; //$$##: number of variance parameters that are time varying.
   else if (dsgemodel_ps->indx_hard_restrict < 200)
      finalpars_ps->n_tvvarpars = 9; //$$##: number of variance parameters that are time varying.
   else if (dsgemodel_ps->indx_hard_restrict < 300)
      finalpars_ps->n_tvvarpars = 8; //$$##: number of variance parameters that are time varying.
   //+
   if (dw_DimA(n_st_var) != finalpars_ps->n_tvvarpars)
      fn_DisplayError(".../*_comfuns.c/CreateTSfinalpars(): number of tv varinaces does not match");
   //+
   finalpars_ps->gsigmaglambdaa_dv = CreateConstantVector_lf(n_st_var[0], 0.0);
   finalpars_ps->gsigmaglambdaz_dv = CreateConstantVector_lf(n_st_var[1], 0.0);
   finalpars_ps->gsigmagnuz_dv     = CreateConstantVector_lf(n_st_var[2], 0.0);
   finalpars_ps->gsigmaglambdaq_dv = CreateConstantVector_lf(n_st_var[3], 0.0);
   finalpars_ps->gsigmagnuq_dv     = CreateConstantVector_lf(n_st_var[4], 0.0);
   finalpars_ps->gsigmagvarphi_dv  = CreateConstantVector_lf(n_st_var[5], 0.0);
   finalpars_ps->gsigmagpsi_dv     = CreateConstantVector_lf(n_st_var[6], 0.0);
   finalpars_ps->gsigmagtheta_dv   = CreateConstantVector_lf(n_st_var[7], 0.0);
   if ((dsgemodel_ps->indx_hard_restrict >= 100) && (dsgemodel_ps->indx_hard_restrict < 200))
      finalpars_ps->gsigmagdelta_dv   = CreateConstantVector_lf(n_st_var[8], 0.0);
   else
      finalpars_ps->gsigmagdelta_dv = NULL;
                                 //$$## hard-coded for [7] or [8].

   return (finalpars_ps);
}
//+
struct TSfinalpars_tag *DestroyTSfinalpars(struct TSfinalpars_tag *finalpars_ps)
{
   if (finalpars_ps)
   {
      //=== The order matters!
      DestroyVector_lf(finalpars_ps->grhoglambdaa_dv);
      DestroyVector_lf(finalpars_ps->grhoglambdaz_dv);
      DestroyVector_lf(finalpars_ps->grhognuz_dv);
      DestroyVector_lf(finalpars_ps->grhoglambdaq_dv);
      DestroyVector_lf(finalpars_ps->grhognuq_dv);
      DestroyVector_lf(finalpars_ps->grhogvarphi_dv);
      DestroyVector_lf(finalpars_ps->grhogpsi_dv);
      DestroyVector_lf(finalpars_ps->grhogtheta_dv);
      DestroyVector_lf(finalpars_ps->grhogdelta_dv);
      DestroyMatrix_lf(finalpars_ps->Qc_dm);
      //+
      DestroyVector_lf(finalpars_ps->gsigmaglambdaa_dv);
      DestroyVector_lf(finalpars_ps->gsigmaglambdaz_dv);
      DestroyVector_lf(finalpars_ps->gsigmagnuz_dv    );
      DestroyVector_lf(finalpars_ps->gsigmaglambdaq_dv);
      DestroyVector_lf(finalpars_ps->gsigmagnuq_dv    );
      DestroyVector_lf(finalpars_ps->gsigmagvarphi_dv );
      DestroyVector_lf(finalpars_ps->gsigmagpsi_dv    );
      DestroyVector_lf(finalpars_ps->gsigmagtheta_dv  );
      DestroyVector_lf(finalpars_ps->gsigmagdelta_dv  );
      //--- To be destroyed last!
      tzDestroy(finalpars_ps);

      return ((struct TSfinalpars_tag *)NULL);
   }
   else  return (finalpars_ps);
}


//------------------------------------
// Converting free parameters (bar transition matrix elements) to final parameters.
//------------------------------------
void Convertphifreepars2finalpars(struct TSdsgemodel_tag *dsgemodel_ps)
{
   struct TSfinalpars_tag *finalpars_ps = dsgemodel_ps->finalpars_ps;
   TSdvector *xphi_dv = dsgemodel_ps->xphi_dv;  //dsgemodel_ps->xphi_dv was refreshed with minpack before this function is called.
   //---
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   int *cum_total_coef = DWRegimeControlVariables_ps->cum_total_coef;
   int *n_st_var = DWRegimeControlVariables_ps->n_st_var;
   int *cum_total_var = DWRegimeControlVariables_ps->cum_total_var;


   //=== $$##: may need to change indx_hard_restrict for the following if loop.
   if (dsgemodel_ps->indx_hard_restrict==0)
   {
      //--- Deep parameters.
      finalpars_ps->gbeta       = 1.0/( 1.0 + 0.01*xphi_dv->v[0] );  //discount factor in the consumer's problem: g_beta = 1.0/(1.0 + 0.01*x_g_beta).
      finalpars_ps->glambdaq    = 1.0 + 0.01*xphi_dv->v[1];  //investment-specific technology growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->glambdastar = 1.0 + 0.01*xphi_dv->v[2];  //consumption growth rate: g_lambda_star = 1.0 + 0.01*x_g_lambda_star.
      finalpars_ps->ggammah     = xphi_dv->v[3];
      finalpars_ps->ggammae     = xphi_dv->v[4];
      finalpars_ps->geta        = xphi_dv->v[5];
      finalpars_ps->gOmega      = xphi_dv->v[6];
   }
   else if (dsgemodel_ps->indx_hard_restrict==1)
   {
      //--- Deep parameters.
      finalpars_ps->gbeta       = 1.0/( 1.0 + 0.01*xphi_dv->v[0] );  //discount factor in the consumer's problem: g_beta = 1.0/(1.0 + 0.01*x_g_beta).
      finalpars_ps->glambdaq    = 1.0 + 0.01*xphi_dv->v[1];  //neutral tech growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->ggammah     = xphi_dv->v[2];
      finalpars_ps->ggammae     = xphi_dv->v[3];
      finalpars_ps->gOmega      = xphi_dv->v[4];
   }
   else if (dsgemodel_ps->indx_hard_restrict==2)
   {
      //--- Deep parameters.
      finalpars_ps->glambdaq    = 1.0 + 0.01*xphi_dv->v[0];  //investment-specific technology growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->glambdastar    = 1.0 + 0.01*xphi_dv->v[1];  //neutral tech growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->ggammah     = xphi_dv->v[2];
      finalpars_ps->ggammae     = xphi_dv->v[3];
      finalpars_ps->geta        = xphi_dv->v[4];
      finalpars_ps->gOmega      = xphi_dv->v[5];
   }
   else if ((dsgemodel_ps->indx_hard_restrict==3) || (dsgemodel_ps->indx_hard_restrict==103) || (dsgemodel_ps->indx_hard_restrict==203))
   {
      //--- Deep parameters.
      finalpars_ps->glambdaq    = 1.0 + 0.01*xphi_dv->v[0];  //investment-specific technology growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->glambdastar    = 1.0 + 0.01*xphi_dv->v[1];  //neutral tech growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->ggammah     = xphi_dv->v[2];
      finalpars_ps->ggammae     = xphi_dv->v[3];
      finalpars_ps->gOmega      = xphi_dv->v[4];
   }
   else if (dsgemodel_ps->indx_hard_restrict==207)
   {
      //--- Deep parameters.
      finalpars_ps->ggammah     = xphi_dv->v[0];
      finalpars_ps->ggammae     = xphi_dv->v[1];
      finalpars_ps->gOmega      = xphi_dv->v[2];
   }
   else if (dsgemodel_ps->indx_hard_restrict==4)
   {
      //--- Deep parameters.
      finalpars_ps->glambdaq    = 1.0 + 0.01*xphi_dv->v[0];  //investment-specific technology growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->ggammah     = xphi_dv->v[1];
      finalpars_ps->ggammae     = xphi_dv->v[2];
      finalpars_ps->gOmega      = xphi_dv->v[3];
   }
   else if ((dsgemodel_ps->indx_hard_restrict==5) || (dsgemodel_ps->indx_hard_restrict==105) || (dsgemodel_ps->indx_hard_restrict==205))
   {
      //--- Deep parameters.
      finalpars_ps->glambdaq    = 1.0 + 0.01*xphi_dv->v[0];  //investment-specific technology growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->glambdastar    = 1.0 + 0.01*xphi_dv->v[1];  //neutral tech growth rate: g_lambda_q = 1.0 + 0.01*x_g_lambda_q.
      finalpars_ps->ggammah     = xphi_dv->v[2];
      finalpars_ps->ggammae     = xphi_dv->v[3];
      finalpars_ps->gOmega      = xphi_dv->v[4];
      finalpars_ps->R           = 1.0 + 0.01*xphi_dv->v[5];
   }
   else  fn_DisplayError(".../*_comfuns.c/Convertphifreepars2finalpars():  indx_hard_restrict in datainp_common.prn is NOT properly specified");


   //------ Derived parameter values.
   #include CInputDerivedSSPars_TAG


   //--- Time-varying coefficients.
   memcpy(finalpars_ps->grhoglambdaa_dv->v, xphi_dv->v+cum_total_coef[0], n_st_coef[0]*sizeof(double));
   memcpy(finalpars_ps->grhoglambdaz_dv->v, xphi_dv->v+cum_total_coef[1], n_st_coef[1]*sizeof(double));
   memcpy(finalpars_ps->grhognuz_dv->v,     xphi_dv->v+cum_total_coef[2], n_st_coef[2]*sizeof(double));
   memcpy(finalpars_ps->grhoglambdaq_dv->v, xphi_dv->v+cum_total_coef[3], n_st_coef[3]*sizeof(double));
   memcpy(finalpars_ps->grhognuq_dv->v,     xphi_dv->v+cum_total_coef[4], n_st_coef[4]*sizeof(double));
   memcpy(finalpars_ps->grhogvarphi_dv->v,  xphi_dv->v+cum_total_coef[5], n_st_coef[5]*sizeof(double));
   memcpy(finalpars_ps->grhogpsi_dv->v,     xphi_dv->v+cum_total_coef[6], n_st_coef[6]*sizeof(double));
   memcpy(finalpars_ps->grhogtheta_dv->v,   xphi_dv->v+cum_total_coef[7], n_st_coef[7]*sizeof(double));
   if ((dsgemodel_ps->indx_hard_restrict >= 100) && (dsgemodel_ps->indx_hard_restrict < 200))
      memcpy(finalpars_ps->grhogdelta_dv->v,   xphi_dv->v+cum_total_coef[8], n_st_coef[8]*sizeof(double));
                   //$$## [7] or [8] must match dw_DimA(n_st_coef) or finalpars_ps->n_tvcoefpars in CreateTSfinalpars().

   //--- Time-varying shock standard deviations.
   memcpy(finalpars_ps->gsigmaglambdaa_dv->v, xphi_dv->v+cum_total_var[0], n_st_var[0]*sizeof(double));
   memcpy(finalpars_ps->gsigmaglambdaz_dv->v, xphi_dv->v+cum_total_var[1], n_st_var[1]*sizeof(double));
   memcpy(finalpars_ps->gsigmagnuz_dv    ->v, xphi_dv->v+cum_total_var[2], n_st_var[2]*sizeof(double));
   memcpy(finalpars_ps->gsigmaglambdaq_dv->v, xphi_dv->v+cum_total_var[3], n_st_var[3]*sizeof(double));
   memcpy(finalpars_ps->gsigmagnuq_dv    ->v, xphi_dv->v+cum_total_var[4], n_st_var[4]*sizeof(double));
   memcpy(finalpars_ps->gsigmagvarphi_dv ->v, xphi_dv->v+cum_total_var[5], n_st_var[5]*sizeof(double));
   memcpy(finalpars_ps->gsigmagpsi_dv    ->v, xphi_dv->v+cum_total_var[6], n_st_var[6]*sizeof(double));
   memcpy(finalpars_ps->gsigmagtheta_dv  ->v, xphi_dv->v+cum_total_var[7], n_st_var[7]*sizeof(double));
   if ((dsgemodel_ps->indx_hard_restrict >= 100) && (dsgemodel_ps->indx_hard_restrict < 200))
      memcpy(finalpars_ps->gsigmagdelta_dv  ->v, xphi_dv->v+cum_total_var[8], n_st_var[8]*sizeof(double));
                   //$$## [7] or [8] must match dw_DimA(n_st_var) or finalpars_ps->n_tvvarpars in CreateTSfinalpars().
}
//------------------------------------
// Checking whether the parameters satisfy the a priori restrictions.
//------------------------------------
int ViolateAPrioriRestrictions(struct TSdsgemodel_tag *dsgemodel_ps)
{
   //Return 1: the parameters violate the a priori restrictions.
   //       0: the parameters satisfy the a priori restrictions.

   int si, sj;
   int indx_violate_apriori_restrictions = 0;  //1: the parameters violate the a priori restrictions.
                                               //0: the parameters satisfy the a priori restrictions.
   TSdvector *xphi_dv = dsgemodel_ps->xphi_dv;
   struct TSfinalpars_tag *finalpars_ps = dsgemodel_ps->finalpars_ps;
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   int *cum_total_coef = DWRegimeControlVariables_ps->cum_total_coef;
   int *n_st_var = DWRegimeControlVariables_ps->n_st_var;
   int *cum_total_var = DWRegimeControlVariables_ps->cum_total_var;



   if (dsgemodel_ps->indx_hard_restrict==0)
   {
      //====== Deep parameters. ======
      if (xphi_dv->v[0] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //discount factor
      if (xphi_dv->v[1] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment tech
      if (xphi_dv->v[2] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //consumption tech
      if (xphi_dv->v[3] < 0.0 || xphi_dv->v[3] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[4] < 0.0 || xphi_dv->v[4] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[5] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //inverse Frisch
      if (xphi_dv->v[6] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment cost
   }
   else if (dsgemodel_ps->indx_hard_restrict==1)
   {
      //====== Deep parameters. ======
      if (xphi_dv->v[0] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //discount factor
      if (xphi_dv->v[1] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment tech
      if (xphi_dv->v[2] < 0.0 || xphi_dv->v[2] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[3] < 0.0 || xphi_dv->v[3] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[4] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment cost
   }
   else if (dsgemodel_ps->indx_hard_restrict==2)
   {
      //====== Deep parameters. ======
      if (xphi_dv->v[0] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment tech
      if (xphi_dv->v[1] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //consumption tech
      if (xphi_dv->v[2] < 0.0 || xphi_dv->v[2] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[3] < 0.0 || xphi_dv->v[3] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[4] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //inverse Frisch
      if (xphi_dv->v[5] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment cost
   }
   else if ((dsgemodel_ps->indx_hard_restrict==3) || (dsgemodel_ps->indx_hard_restrict==103) || (dsgemodel_ps->indx_hard_restrict==203))
   {
      //====== Deep parameters. ======
      if (xphi_dv->v[0] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment tech
      if (xphi_dv->v[1] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //consumption tech
      if (xphi_dv->v[2] < 0.0 || xphi_dv->v[2] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[3] < 0.0 || xphi_dv->v[3] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[4] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment cost
   }
   else if (dsgemodel_ps->indx_hard_restrict==207)
   {
      //====== Deep parameters. ======
      if (xphi_dv->v[0] < 0.0 || xphi_dv->v[2] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[1] < 0.0 || xphi_dv->v[3] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[2] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment cost
   }
   else if (dsgemodel_ps->indx_hard_restrict==4)
   {
      //====== Deep parameters. ======
      if (xphi_dv->v[0] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment tech
      if (xphi_dv->v[1] < 0.0 || xphi_dv->v[1] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[2] < 0.0 || xphi_dv->v[2] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[3] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment cost
   }
   else if ((dsgemodel_ps->indx_hard_restrict==5) || (dsgemodel_ps->indx_hard_restrict==105) || (dsgemodel_ps->indx_hard_restrict==205))
   {
      //====== Deep parameters. ======
      if (xphi_dv->v[0] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment tech
      if (xphi_dv->v[1] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //consumption tech
      if (xphi_dv->v[2] < 0.0 || xphi_dv->v[2] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[3] < 0.0 || xphi_dv->v[3] >= 1.0)  return (indx_violate_apriori_restrictions = 1);  //habit
      if (xphi_dv->v[4] < 0.0)  return (indx_violate_apriori_restrictions = 1);  //investment cost
      if (xphi_dv->v[5] <= 0.0)  return (indx_violate_apriori_restrictions = 1);  //risk-free loan rate
   }


   //--- TV coefficients
   for (sj=finalpars_ps->n_tvcoefpars-1; sj>=0; sj--)
      for (si=n_st_coef[sj]-1; si>=0; si--)
         if (xphi_dv->v[cum_total_coef[sj]+si] < 0.0)  return (indx_violate_apriori_restrictions = 1);

   //--- Shock variances.
   for (sj=finalpars_ps->n_tvvarpars-1; sj>=0; sj--)
      for (si=n_st_var[sj]-1; si>=0; si--)
         if (xphi_dv->v[cum_total_var[sj]+si] <= 0.0)  return (indx_violate_apriori_restrictions = 1);   //standard deviation.

   //----------------------------------------------------
   //- Restrictions that must be satisfied.
   //- With qlLeOy, yk, ylLhOy calibrated, these restrictions should be satisfied automatically.
   //----------------------------------------------------
   if (finalpars_ps->glambdastar >= (finalpars_ps->R - SQRTEPSILON))  //This restriction means that it must hold for gbeta*(1+lambda_a) < 1.
      return (indx_violate_apriori_restrictions = 1);
   if (finalpars_ps->gtheta >= ((1.0-finalpars_ps->gbeta)/(finalpars_ps->gbeta*finalpars_ps->glambdaa) - SQRTEPSILON))
      return (indx_violate_apriori_restrictions = 1);
   //--- We should allow glambdaz to be less than 1.0, I think.  12 June 2009.
   //if (finalpars_ps->glambdaz < 1.0 -  SQRTEPSILON)
   //   return (indx_violate_apriori_restrictions = 1);
   if (finalpars_ps->glambdaz < 0.0 +  SQRTEPSILON)
      return (indx_violate_apriori_restrictions = 1);



   return (indx_violate_apriori_restrictions);
}


//------------------------------------
// Creates prior hyperparameter values for the prior of model-specific parameters (bar transition matrix elements)
//    for both the constant-parameter and the regime-switching models.
//------------------------------------
struct TSpriorhyperpars_tag *CreateTSpriorhyperpars()
{
   //The prior hyperparameters are MANUALLY changed or keyed in.
   //
   //See Zhe Notes.

   struct TSpriorhyperpars_tag *priorhyperpars_ps = tzMalloc(1, struct TSpriorhyperpars_tag);


   //====== Deep parameters. ======
   priorhyperpars_ps->a_x_gbeta = 1.5832;
   priorhyperpars_ps->b_x_gbeta = 1.0126;
   //+
   priorhyperpars_ps->a_x_glambdaq = 1.8611;
   priorhyperpars_ps->b_x_glambdaq = 3.0112;
   priorhyperpars_ps->a_x_glambdastar = 1.8611;
   priorhyperpars_ps->b_x_glambdastar = 3.0112;
   //+
   priorhyperpars_ps->a_ggammah = 1.0;
   priorhyperpars_ps->b_ggammah = 2.0;
   priorhyperpars_ps->a_ggammae = 1.0;
   priorhyperpars_ps->b_ggammae = 2.0;
   //+
   priorhyperpars_ps->a_geta = 1.0;
   priorhyperpars_ps->b_geta = 1.0;
   priorhyperpars_ps->a_gOmega = 1.0;
   priorhyperpars_ps->b_gOmega = 0.5;
   //
   priorhyperpars_ps->a_x_R = 0.7371;
   priorhyperpars_ps->b_x_R = 0.3078;



   //=== Exogenous processes
   //--- AR persistence
   // priorhyperpars_ps->a_grhoglambdaa = 1.0;
   // priorhyperpars_ps->a_grhoglambdaz = 1.0;
   // priorhyperpars_ps->a_grhognuz = 1.0;
   // priorhyperpars_ps->a_grhoglambdaq = 1.0;
   // priorhyperpars_ps->a_grhognuq = 1.0;
   // priorhyperpars_ps->a_grhogvarphi = 1.0;
   // priorhyperpars_ps->a_grhogpsi = 1.0;
   // priorhyperpars_ps->a_grhogtheta = 1.0;
   // //+
   // priorhyperpars_ps->b_grhoglambdaa = 2.0;
   // priorhyperpars_ps->b_grhoglambdaz = 2.0;
   // priorhyperpars_ps->b_grhognuz = 2.0;
   // priorhyperpars_ps->b_grhoglambdaq = 2.0;
   // priorhyperpars_ps->b_grhognuq = 2.0;
   // priorhyperpars_ps->b_grhogvarphi = 2.0;
   // priorhyperpars_ps->b_grhogpsi = 2.0;
   // priorhyperpars_ps->b_grhogtheta = 2.0;
   //+
   priorhyperpars_ps->a_rho_all = 1.0;
   priorhyperpars_ps->b_rho_all = 2.0;

   //--- Standard deviations of all shocks.
   // priorhyperpars_ps->a_sigma_all = 0.5675;    //Inverse-Gamma prior
   // priorhyperpars_ps->b_sigma_all = 0.0021;    //Inverse-Gamma prior
   //   priorhyperpars_ps->a_sigma_all = 0.4979;    //Inverse-Gamma prior
   //   priorhyperpars_ps->b_sigma_all = 0.0002;    //Inverse-Gamma prior
   //   priorhyperpars_ps->a_sigma_all = 0.5943;    //Inverse-Gamma prior
   //   priorhyperpars_ps->b_sigma_all = 0.0011;    //Inverse-Gamma prior
   //// priorhyperpars_ps->a_sigma_all = 3.5428e-01;    //Inverse-Gamma prior [0.0001 1.0]
   //// priorhyperpars_ps->b_sigma_all = 1.5342e-04;    //Inverse-Gamma prior [0.0001 1.0]
   priorhyperpars_ps->a_sigma_all = 3.260961642159758e-01;    //Inverse-Gamma prior [0.0001 2.0]
   priorhyperpars_ps->b_sigma_all = 1.451789062500018e-04;    //Inverse-Gamma prior [0.0001 2.0]
     
   return (priorhyperpars_ps);
}
//---
//struct TSpriorhyperpars_tag *DestroyTSpriorhyperpars(struct TSpriorhyperpars_tag *priorhyperpars_ps)
//{
//   if (priorhyperpars_ps)
//   {
//      //+  To be destroyed last!
//      tzDestroy(priorhyperpars_ps);
//
//      return ((struct TSpriorhyperpars_tag *)NULL);
//   }
//   else  return (priorhyperpars_ps);
//}


//------------------------------------
// log prior density for both the constant-parameter and regime-switching models.
//------------------------------------

double logpriordensity(struct TStateModel_tag *smodel_ps)
{
   //The prior parameters are MANUALLY changed or keyed in.
   //
   //See pp. 49-50 in LWZ Model II NOTES.

   int si, sj;
   double logprior;
   //
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   int *cum_total_coef = DWRegimeControlVariables_ps->cum_total_coef;
   int *n_st_var = DWRegimeControlVariables_ps->n_st_var;
   int *cum_total_var = DWRegimeControlVariables_ps->cum_total_var;
   //
   TSdvector *xphi_dv = dsgemodel_ps->xphi_dv;
   struct TSfinalpars_tag *finalpars_ps = dsgemodel_ps->finalpars_ps;
   struct TSpriorhyperpars_tag *priorhyperpars_ps = dsgemodel_ps->priorhyperpars_ps;


   //--- The following calls are absolutely NECESSARY becasue (1) the prior may be called first before
   //---   logTimetCondLH() is called, and after the prior is called, dsgemodel_ps->refreshed will be set
   //---   to 1, so unless we call RefreshEverything(), gensys will not be called (tragically)
   //---   (2) RefreshEverthing() will be called only ONCE for the integrated-out-over-regimes LH from tbase0=0 to tbase0=T-1.
   //---   until the DW's reset function ThetaChanged() is called, which calls tz_thetaChanged() and tz_TransitionMatrixChanged().
   //---   DW's code will call tz_thetaChanged() and tz_TransitionMatrixChanged() only when
   //---   parameters are changed (e.g., when inpt=T (reaches the end of the sample)
   //---   and new parameters are given and the forward recursion starts from the beginning again.).
   if (!dsgemodel_ps->refreshed)
   {
      if (!RefreshEverything(smodel_ps))  return (-NEARINFINITY);
   }
   else
   {
      if (!dsgemodel_ps->ValidParameters || !dsgemodel_ps->ValidGensysSolution || !dsgemodel_ps->ValidInitial_z10_P10)  return (-NEARINFINITY);
   }



   logprior = 0.0;

   if (dsgemodel_ps->indx_hard_restrict==0)
   {
      logprior += tz_loggammapdf(xphi_dv->v[0], priorhyperpars_ps->a_x_gbeta, priorhyperpars_ps->b_x_gbeta);
      logprior += tz_loggammapdf(xphi_dv->v[1], priorhyperpars_ps->a_x_glambdaq, priorhyperpars_ps->b_x_glambdaq);
      logprior += tz_loggammapdf(xphi_dv->v[2], priorhyperpars_ps->a_x_glambdastar, priorhyperpars_ps->b_x_glambdastar);
      logprior += tz_logbetapdf(xphi_dv->v[3], priorhyperpars_ps->a_ggammah, priorhyperpars_ps->b_ggammah);
      logprior += tz_logbetapdf(xphi_dv->v[4], priorhyperpars_ps->a_ggammae, priorhyperpars_ps->b_ggammae);
      logprior += tz_loggammapdf(xphi_dv->v[5], priorhyperpars_ps->a_geta, priorhyperpars_ps->b_geta);
      logprior += tz_loggammapdf(xphi_dv->v[6], priorhyperpars_ps->a_gOmega, priorhyperpars_ps->b_gOmega);
   }
   else if (dsgemodel_ps->indx_hard_restrict==1)
   {
      //====== Deep parameters. ======
      logprior += tz_loggammapdf(xphi_dv->v[0], priorhyperpars_ps->a_x_gbeta, priorhyperpars_ps->b_x_gbeta);
      logprior += tz_loggammapdf(xphi_dv->v[1], priorhyperpars_ps->a_x_glambdaq, priorhyperpars_ps->b_x_glambdaq);
      logprior += tz_logbetapdf(xphi_dv->v[2], priorhyperpars_ps->a_ggammah, priorhyperpars_ps->b_ggammah);
      logprior += tz_logbetapdf(xphi_dv->v[3], priorhyperpars_ps->a_ggammae, priorhyperpars_ps->b_ggammae);
      logprior += tz_loggammapdf(xphi_dv->v[4], priorhyperpars_ps->a_gOmega, priorhyperpars_ps->b_gOmega);
   }
   else if (dsgemodel_ps->indx_hard_restrict==2)
   {
      logprior += tz_loggammapdf(xphi_dv->v[0], priorhyperpars_ps->a_x_glambdaq, priorhyperpars_ps->b_x_glambdaq);
      logprior += tz_loggammapdf(xphi_dv->v[1], priorhyperpars_ps->a_x_glambdastar, priorhyperpars_ps->b_x_glambdastar);
      logprior += tz_logbetapdf(xphi_dv->v[2], priorhyperpars_ps->a_ggammah, priorhyperpars_ps->b_ggammah);
      logprior += tz_logbetapdf(xphi_dv->v[3], priorhyperpars_ps->a_ggammae, priorhyperpars_ps->b_ggammae);
      logprior += tz_loggammapdf(xphi_dv->v[4], priorhyperpars_ps->a_geta, priorhyperpars_ps->b_geta);
      logprior += tz_loggammapdf(xphi_dv->v[5], priorhyperpars_ps->a_gOmega, priorhyperpars_ps->b_gOmega);
   }
   else if ((dsgemodel_ps->indx_hard_restrict==3) || (dsgemodel_ps->indx_hard_restrict==103) || (dsgemodel_ps->indx_hard_restrict==203))
   {
      logprior += tz_loggammapdf(xphi_dv->v[0], priorhyperpars_ps->a_x_glambdaq, priorhyperpars_ps->b_x_glambdaq);
      logprior += tz_loggammapdf(xphi_dv->v[1], priorhyperpars_ps->a_x_glambdastar, priorhyperpars_ps->b_x_glambdastar);
      logprior += tz_logbetapdf(xphi_dv->v[2], priorhyperpars_ps->a_ggammah, priorhyperpars_ps->b_ggammah);
      logprior += tz_logbetapdf(xphi_dv->v[3], priorhyperpars_ps->a_ggammae, priorhyperpars_ps->b_ggammae);
      logprior += tz_loggammapdf(xphi_dv->v[4], priorhyperpars_ps->a_gOmega, priorhyperpars_ps->b_gOmega);
   }
   else if (dsgemodel_ps->indx_hard_restrict==207)
   {
      logprior += tz_logbetapdf(xphi_dv->v[0], priorhyperpars_ps->a_ggammah, priorhyperpars_ps->b_ggammah);
      logprior += tz_logbetapdf(xphi_dv->v[1], priorhyperpars_ps->a_ggammae, priorhyperpars_ps->b_ggammae);
      logprior += tz_loggammapdf(xphi_dv->v[2], priorhyperpars_ps->a_gOmega, priorhyperpars_ps->b_gOmega);
   }
   else if (dsgemodel_ps->indx_hard_restrict==4)
   {
      logprior += tz_loggammapdf(xphi_dv->v[0], priorhyperpars_ps->a_x_glambdaq, priorhyperpars_ps->b_x_glambdaq);
      logprior += tz_logbetapdf(xphi_dv->v[1], priorhyperpars_ps->a_ggammah, priorhyperpars_ps->b_ggammah);
      logprior += tz_logbetapdf(xphi_dv->v[2], priorhyperpars_ps->a_ggammae, priorhyperpars_ps->b_ggammae);
      logprior += tz_loggammapdf(xphi_dv->v[3], priorhyperpars_ps->a_gOmega, priorhyperpars_ps->b_gOmega);
   }
   else if ((dsgemodel_ps->indx_hard_restrict==5) || (dsgemodel_ps->indx_hard_restrict==105) || (dsgemodel_ps->indx_hard_restrict==205))
   {
      logprior += tz_loggammapdf(xphi_dv->v[0], priorhyperpars_ps->a_x_glambdaq, priorhyperpars_ps->b_x_glambdaq);
      logprior += tz_loggammapdf(xphi_dv->v[1], priorhyperpars_ps->a_x_glambdastar, priorhyperpars_ps->b_x_glambdastar);
      logprior += tz_logbetapdf(xphi_dv->v[2], priorhyperpars_ps->a_ggammah, priorhyperpars_ps->b_ggammah);
      logprior += tz_logbetapdf(xphi_dv->v[3], priorhyperpars_ps->a_ggammae, priorhyperpars_ps->b_ggammae);
      logprior += tz_loggammapdf(xphi_dv->v[4], priorhyperpars_ps->a_gOmega, priorhyperpars_ps->b_gOmega);
      logprior += tz_loggammapdf(xphi_dv->v[5], priorhyperpars_ps->a_x_R, priorhyperpars_ps->b_x_R);
   }


   //--- TV coefficients.
   for (sj=finalpars_ps->n_tvvarpars-1; sj>=0; sj--)
      for (si=n_st_coef[sj]-1; si>=0; si--)
         logprior += tz_logbetapdf(xphi_dv->v[cum_total_coef[sj]+si], priorhyperpars_ps->a_rho_all, priorhyperpars_ps->b_rho_all);

   //--- Shock standard deviations.
   for (sj=finalpars_ps->n_tvvarpars-1; sj>=0; sj--)
      for (si=n_st_var[sj]-1; si>=0; si--)
         logprior += tz_loginversegammapdf(xphi_dv->v[cum_total_var[sj]+si], priorhyperpars_ps->a_sigma_all, priorhyperpars_ps->b_sigma_all); //standard deviations of different shocks looping through.

   return (logprior);
}


void Draw_xphi_dvFromPrior_Specific(struct TStateModel_tag *smodel_ps)
{
   int si, sj;
   //--- Taken out of Waggoner's model.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   int *cum_total_coef = DWRegimeControlVariables_ps->cum_total_coef;
   int *n_st_var = DWRegimeControlVariables_ps->n_st_var;
   int *cum_total_var = DWRegimeControlVariables_ps->cum_total_var;
   //
   TSdvector *xphi_dv = dsgemodel_ps->xphi_dv;
   struct TSfinalpars_tag *finalpars_ps = dsgemodel_ps->finalpars_ps;
   struct TSpriorhyperpars_tag *priorhyperpars_ps = dsgemodel_ps->priorhyperpars_ps;

   if (dsgemodel_ps->indx_hard_restrict==0)
   {
      xphi_dv->v[0]  = GammaDouble( priorhyperpars_ps->a_x_gbeta      , priorhyperpars_ps->b_x_gbeta      );
      xphi_dv->v[1]  = GammaDouble( priorhyperpars_ps->a_x_glambdaq   , priorhyperpars_ps->b_x_glambdaq   );
      xphi_dv->v[2]  = GammaDouble( priorhyperpars_ps->a_x_glambdastar, priorhyperpars_ps->b_x_glambdastar);
      xphi_dv->v[3]  = betarnd(     priorhyperpars_ps->a_ggammah      , priorhyperpars_ps->b_ggammah      );
      xphi_dv->v[4]  = betarnd(     priorhyperpars_ps->a_ggammae      , priorhyperpars_ps->b_ggammae      );
      xphi_dv->v[5]  = GammaDouble( priorhyperpars_ps->a_geta         , priorhyperpars_ps->b_geta         );
      xphi_dv->v[6]  = GammaDouble( priorhyperpars_ps->a_gOmega       , priorhyperpars_ps->b_gOmega       );
   }
   else if (dsgemodel_ps->indx_hard_restrict==1)
   {
      //====== Deep parameters. ======
      xphi_dv->v[0]  = GammaDouble( priorhyperpars_ps->a_x_gbeta   , priorhyperpars_ps->b_x_gbeta   );
      xphi_dv->v[1]  = GammaDouble( priorhyperpars_ps->a_x_glambdaq, priorhyperpars_ps->b_x_glambdaq);
      xphi_dv->v[2]  = betarnd(     priorhyperpars_ps->a_ggammah   , priorhyperpars_ps->b_ggammah   );
      xphi_dv->v[3]  = betarnd(     priorhyperpars_ps->a_ggammae   , priorhyperpars_ps->b_ggammae   );
      xphi_dv->v[4]  = GammaDouble( priorhyperpars_ps->a_gOmega    , priorhyperpars_ps->b_gOmega    );
   }
   else if (dsgemodel_ps->indx_hard_restrict==2)
   {
      xphi_dv->v[0]  = GammaDouble( priorhyperpars_ps->a_x_glambdaq   , priorhyperpars_ps->b_x_glambdaq   );
      xphi_dv->v[1]  = GammaDouble( priorhyperpars_ps->a_x_glambdastar, priorhyperpars_ps->b_x_glambdastar);
      xphi_dv->v[2]  = betarnd(     priorhyperpars_ps->a_ggammah      , priorhyperpars_ps->b_ggammah      );
      xphi_dv->v[3]  = betarnd(     priorhyperpars_ps->a_ggammae      , priorhyperpars_ps->b_ggammae      );
      xphi_dv->v[4]  = GammaDouble( priorhyperpars_ps->a_geta         , priorhyperpars_ps->b_geta         );
      xphi_dv->v[5]  = GammaDouble( priorhyperpars_ps->a_gOmega     , priorhyperpars_ps->b_gOmega       );
   }
   else if ((dsgemodel_ps->indx_hard_restrict==3) || (dsgemodel_ps->indx_hard_restrict==103) || (dsgemodel_ps->indx_hard_restrict==203))
   {
      xphi_dv->v[0]  = GammaDouble( priorhyperpars_ps->a_x_glambdaq   , priorhyperpars_ps->b_x_glambdaq   );
      xphi_dv->v[1]  = GammaDouble( priorhyperpars_ps->a_x_glambdastar, priorhyperpars_ps->b_x_glambdastar);
      xphi_dv->v[2]  = betarnd(     priorhyperpars_ps->a_ggammah      , priorhyperpars_ps->b_ggammah      );
      xphi_dv->v[3]  = betarnd(     priorhyperpars_ps->a_ggammae      , priorhyperpars_ps->b_ggammae      );
      xphi_dv->v[4]  = GammaDouble( priorhyperpars_ps->a_gOmega       , priorhyperpars_ps->b_gOmega       );
   }
   else if (dsgemodel_ps->indx_hard_restrict==207)
   {
      xphi_dv->v[0]  = betarnd(     priorhyperpars_ps->a_ggammah      , priorhyperpars_ps->b_ggammah      );
      xphi_dv->v[1]  = betarnd(     priorhyperpars_ps->a_ggammae      , priorhyperpars_ps->b_ggammae      );
      xphi_dv->v[2]  = GammaDouble( priorhyperpars_ps->a_gOmega       , priorhyperpars_ps->b_gOmega       );
   }
   else if (dsgemodel_ps->indx_hard_restrict==4)
   {
      xphi_dv->v[0]  = GammaDouble( priorhyperpars_ps->a_x_glambdaq, priorhyperpars_ps->b_x_glambdaq);
      xphi_dv->v[1]  = betarnd(     priorhyperpars_ps->a_ggammah   , priorhyperpars_ps->b_ggammah   );
      xphi_dv->v[2]  = betarnd(     priorhyperpars_ps->a_ggammae   , priorhyperpars_ps->b_ggammae   );
      xphi_dv->v[3]  = GammaDouble( priorhyperpars_ps->a_gOmega    , priorhyperpars_ps->b_gOmega    );
   }
   else if ((dsgemodel_ps->indx_hard_restrict==5) || (dsgemodel_ps->indx_hard_restrict==105) || (dsgemodel_ps->indx_hard_restrict==205))
   {
      xphi_dv->v[0]  = GammaDouble( priorhyperpars_ps->a_x_glambdaq   , priorhyperpars_ps->b_x_glambdaq   );
      xphi_dv->v[1]  = GammaDouble( priorhyperpars_ps->a_x_glambdastar, priorhyperpars_ps->b_x_glambdastar);
      xphi_dv->v[2]  = betarnd(     priorhyperpars_ps->a_ggammah      , priorhyperpars_ps->b_ggammah      );
      xphi_dv->v[3]  = betarnd(     priorhyperpars_ps->a_ggammae      , priorhyperpars_ps->b_ggammae      );
      xphi_dv->v[4]  = GammaDouble( priorhyperpars_ps->a_gOmega       , priorhyperpars_ps->b_gOmega       );
      xphi_dv->v[5]  = GammaDouble( priorhyperpars_ps->a_x_R          , priorhyperpars_ps->b_x_R          );
   }

   //===== The following is all automated. ============
   //--- TV coefficients.
   for (sj=finalpars_ps->n_tvvarpars-1; sj>=0; sj--)
      for (si=n_st_coef[sj]-1; si>=0; si--)
         xphi_dv->v[cum_total_coef[sj]+si] = betarnd(priorhyperpars_ps->a_rho_all, priorhyperpars_ps->b_rho_all);

   //--- Shock standard deviations.
   for (sj=finalpars_ps->n_tvvarpars-1; sj>=0; sj--)
      for (si=n_st_var[sj]-1; si>=0; si--)
         xphi_dv->v[cum_total_var[sj]+si] = InverseGammaDouble(priorhyperpars_ps->a_sigma_all, priorhyperpars_ps->b_sigma_all); //standard deviations of different shocks looping through.
}

// Deal with constraints
//struct TSModelConstraints_tag *modelconstratints_ps = CreateTSModelConstraints(20, 0, 0, &confun);

