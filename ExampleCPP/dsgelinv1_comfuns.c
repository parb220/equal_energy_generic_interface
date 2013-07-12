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

   fprintf(tzGetDebugFile(), "\n-----------\n");
   fprintf(tzGetDebugFile(), "logLkelihood(s_t=%d, inpt=%d):\n", s_t, inpt);
   fprintf(tzGetDebugFile(), " %10.5f\n", jnk);

   fflush(tzGetDebugFile());

   t_loglht = -(0.5*ny)*LOG2PI - 0.5*logdeterminant(Dtdata_dm) - 0.5*VectorDotVector(wny_dv, etdata_dv);
   fprintf(tzGetDebugFile(), " %10.5f\n", t_loglht);

   fprintf(tzGetDebugFile(), "%%st=%d, inpt=%d, and sti=%d\n", st, inpt, sti);

   fprintf(tzGetDebugFile(), "wP0_dv:\n");
   WriteVector(tzGetDebugFile(), wP0_dv, " %10.5f ");
   fprintf(tzGetDebugFile(), "Vt_dc->C[sti_v=%d]:\n", sti_v);
   WriteMatrix(tzGetDebugFile(), Vt_dc->C[sti_v], " %10.5f ");

   fflush(tzGetDebugFile());
}
/**/

#include "dsgelinv1_comfuns.h"


static struct TSDWRegimeControlVariables_tag *DestroyTSDWRegimeControlVariables(struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps);
//static void Convertphifreepars2finalpars(struct TSdsgemodel_tag *dsgemodel_ps);  //Converting free parameters (bar transition matrix elements) to the final (DSGE) meaningful parameters.
static void RunningGensys_allcases(struct TStateModel_tag *smodel_ps);   //MANUALLY changed for different ny, nz, and ngensys imported from the input data file datainp_TAG.prn.
static void Refresh_kalfilms_allcases(struct TStateModel_tag *smodel_ps);
//static int ViolateAPrioriRestrictions(struct TSdsgemodel_tag *dsgemodel_ps);
static struct TSpriorhyperpars_tag *DestroyTSpriorhyperpars(struct TSpriorhyperpars_tag *priorhyperpars_ps);
static double loglh_const(struct TStateModel_tag *smodel_ps);


//************************************************************
//** Creating the specific model structure and other related structures.
//************************************************************
struct TSdsgemodel_tag *CreateTSdsgemodel(FILE *fptr_rawdata, FILE *fptr_common, FILE *fptr_specific, FILE *fptr_markov, int indxStartValuesForMin)
{
   //Outputs:
   //  dsgemodel_ps: model-specific structure with filled arguments.
   //Inputs:
   //  fptr_rawdata:  pointer to raw data input file.
   //  fptr_common:  pointer to common setup input file.
   //  fptr_specific:  pointer to model-specific setup input file.

   int nrows, ncolse;
   int nvec;
   int loc_row_samplebegin;
   int i_;
   //===
   TSdmatrix *Datarawe_dm = NULL;


   //~~~ Creating the structure and initializing the NULL pointers.
   struct TSdsgemodel_tag *dsgemodel_ps = tzMalloc(1, struct TSdsgemodel_tag);
   dsgemodel_ps->modeltag_cv = NULL;
   //
   dsgemodel_ps->priorhyperpars_ps = NULL;
   dsgemodel_ps->DWRegimeControlVariables_ps = NULL;
   dsgemodel_ps->finalpars_ps = NULL;
   //
   dsgemodel_ps->gensys_ps = NULL;
   //
   dsgemodel_ps->imsl_lh_locs_iv = NULL;
   dsgemodel_ps->imsl_lh_vals_dv = NULL;
   dsgemodel_ps->imsl_rh_vals_dv = NULL;
   //+
   dsgemodel_ps->hard7_imsl_lh_locs_iv = NULL;
   dsgemodel_ps->hard7_imsl_lh_vals_dv = NULL;
   dsgemodel_ps->hard7_imsl_rh_vals_dv = NULL;
   //+
   dsgemodel_ps->hard0_simple_locs_iv = NULL;
   dsgemodel_ps->hard0_simple_lowvals_dv = NULL;
   dsgemodel_ps->hard0_simple_highvals_dv = NULL;
   //+
   dsgemodel_ps->hard1_simple_locs_iv = NULL;
   dsgemodel_ps->hard1_simple_lowvals_dv = NULL;
   dsgemodel_ps->hard1_simple_highvals_dv = NULL;
   //+
   dsgemodel_ps->hard2_simple_locs_iv = NULL;
   dsgemodel_ps->hard2_simple_lowvals_dv = NULL;
   dsgemodel_ps->hard2_simple_highvals_dv = NULL;
   //
   dsgemodel_ps->Dataraw_dm = NULL;
   dsgemodel_ps->Datarawtran_dm = NULL;
   dsgemodel_ps->BeginEndData_iv = NULL;
   //+
   dsgemodel_ps->SampleData_dm = NULL;
   dsgemodel_ps->SampleDatatran_dm = NULL;
   dsgemodel_ps->BeginEndSample_iv = NULL;
   dsgemodel_ps->SelectData_iv = NULL;
   dsgemodel_ps->dates_dv = NULL;
   //
   dsgemodel_ps->xphi_dv = NULL;
   dsgemodel_ps->Omega_dm = NULL;
   dsgemodel_ps->Hessian_dm = NULL;
   //--- Output arguemnts.
   dsgemodel_ps->dw_ProbStates_dv = NULL;
   dsgemodel_ps->PostProbS_dm = NULL;
   dsgemodel_ps->BaseProbS_dc = NULL;
   dsgemodel_ps->strshks_dm = NULL;



   //--- Taking key input values.
   dsgemodel_ps->indxStartValuesForMin = indxStartValuesForMin;
   //--- Key initialization.
   dsgemodel_ps->refreshed = 0;  //Meaning: we need to refresh everying at the beginning of the program.
   dsgemodel_ps->modeltag_cv = CreateVector_c(TAGLEN);    //allocated to 512 for a possible long name or tag.


   //--------- Reads data from input data files. The order matters! ---------
   //=== Reads integers.
   //--- Common setup inputs
   if ( !fn_SetFilePosition(fptr_common, "//== indxEstFinder ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->indxEstFinder) != 1 )
      dsgemodel_ps->indxEstFinder = 0;  //Default   fn_DisplayError("Fatal Error: check the below the line //== indxEstFinder ==// in the common input data file datainp_common*.prn");
   if ( !fn_SetFilePosition(fptr_common, "//== indx_hard_restrict ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->indx_hard_restrict) != 1 )
      dsgemodel_ps->indx_hard_restrict = 0;  //Default: Model 0.
   ////if ( !fn_SetFilePosition(fptr_common, "//== indxSimulateMCMC ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->indxSimulateMCMC) != 1 )
   ////   dsgemodel_ps->indxSimulateMCMC = 0; //Default 
   if ( !fn_SetFilePosition(fptr_common, "//== randomseed ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->randomseed) != 1 )
      dsgemodel_ps->randomseed = 0; //Default: random seed according to the clock time.
   if ( !fn_SetFilePosition(fptr_common, "//== nDrawsInitsp ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->nDrawsInitsp) != 1 )
      dsgemodel_ps->nDrawsInitsp = 100; //Default.
   if ( !fn_SetFilePosition(fptr_common, "//== nDrawsFromPrior ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->nDrawsFromPrior) != 1 )
      dsgemodel_ps->nDrawsFromPrior = 50000; //Default.
   if ( !fn_SetFilePosition(fptr_common, "//== flag_ComputingHessian ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->flag_ComputingHessian) != 1 )
      dsgemodel_ps->flag_ComputingHessian = 0; //Default: no Hessian computed.
   if ( !fn_SetFilePosition(fptr_common, "//== flag_ComputingHist ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->flag_ComputingHist) != 1 )
      dsgemodel_ps->flag_ComputingHist = 0; //Default: no historical decompositions computed.
   if ( !fn_SetFilePosition(fptr_common, "//== Use_IMSL_opt ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->Use_IMSL_opt) != 1 )
      dsgemodel_ps->Use_IMSL_opt = 0; //Default: no IMSL optimization use.
   if ( !fn_SetFilePosition(fptr_common, "//== Use_NPSOL_opt ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->Use_NPSOL_opt) != 1 )
      dsgemodel_ps->Use_NPSOL_opt = 0; //Default: no IMSL optimization use.
   //--- Model-specific setup inputs
   if ( !fn_SetFilePosition(fptr_specific, "//== ny ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->ny) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== ny ==// in the model-dependent input data file datainp_TAG.prn");
   if ( !fn_SetFilePosition(fptr_specific, "//== nlags ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->nlags) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== nlags ==// in the model-dependent input data file datainp_TAG.prn");
   if ( !fn_SetFilePosition(fptr_specific, "//== nz ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->nz) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== nz ==// in the model-dependent input data file datainp_TAG.prn");
   if ( !fn_SetFilePosition(fptr_specific, "//== ngensys ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->ngensys) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== ngensys ==// in the model-dependent input data file datainp_TAG.prn");
   if ( !fn_SetFilePosition(fptr_specific, "//== nzbase ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->nzbase) != 1 )
      dsgemodel_ps->nzbase = dsgemodel_ps->ngensys;
   if ( !fn_SetFilePosition(fptr_specific, "//== nu ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->nu) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== nu ==// in the model-dependent input data file datainp_TAG.prn");
   if ( !fn_SetFilePosition(fptr_specific, "//== ne ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->ne) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== ne ==// in the model-dependent input data file datainp_TAG.prn");
   if ( !fn_SetFilePosition(fptr_specific, "//== nshocks ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->nshocks) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== nshocks ==// in the model-dependent input data file datainp_TAG.prn");
   if ( !fn_SetFilePosition(fptr_specific, "//== nshocks_kalman ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->nshocks_kalman) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== nshocks_kalman ==// in the model-dependent input data file datainp_TAG.prn");
   if ( !fn_SetFilePosition(fptr_specific, "//== nexperrs ==//") || fscanf(fptr_specific, " %d ", &dsgemodel_ps->nexperrs) != 1 )
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the below the line //== nexperrs ==// in the model-dependent input data file datainp_TAG.prn");



   if (dsgemodel_ps->Use_IMSL_opt)
   {
      if ( !fn_SetFilePosition(fptr_common, "//== imsl_neqs ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->imsl_neqs) != 1 )
         dsgemodel_ps->imsl_neqs = 0; //Default: no equality restrictions.
      if ( !fn_SetFilePosition(fptr_common, "//== imsl_ncons ==//") || fscanf(fptr_common, " %d ", &dsgemodel_ps->imsl_ncons) != 1 )
         dsgemodel_ps->imsl_ncons = 0;  //Default: no nonlinear restrictions.
   }
  
   //=== Reads integer vectors.
   //--- Raw data inputs
   if (fn_SetFilePosition(fptr_rawdata, "//== BeginEndData_iv ==//"))
   {
      if ( fscanf(fptr_rawdata, " %d ", &nvec) != 1)
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== BeginEndData_iv ==// in the raw data file dataraw_*.prn");
      dsgemodel_ps->BeginEndData_iv = CreateVector_int(nvec);
      if ( !ReadVector_int(fptr_rawdata, dsgemodel_ps->BeginEndData_iv) )
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== BeginEndData_iv ==// in the raw data file dataraw_*.prn");
      dsgemodel_ps->BeginEndData_iv->flag = V_DEF;
   }
   else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== BeginEndData_iv ==// in the raw data file dataraw_*.prn does not exist");
   //--- Sample and variable selections
   if (fn_SetFilePosition(fptr_specific, "//== BeginEndSample_iv ==//"))
   {
      if ( fscanf(fptr_specific, " %d ", &nvec) != 1)
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== BeginEndSample_iv ==// in the model-dependent input data file datainp_TAG.prn");
      dsgemodel_ps->BeginEndSample_iv = CreateVector_int(nvec);
      if ( !ReadVector_int(fptr_specific, dsgemodel_ps->BeginEndSample_iv) )
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data vector after the first row below the line //== BeginEndSample_iv ==// in the model-dependent input data file datainp_TAG.prn");
      dsgemodel_ps->BeginEndSample_iv->flag = V_DEF;
   }
   else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== BeginEndSample_iv ==// in the model-dependent input data file datainp_TAG.prn does not exist");
   //+
   if (fn_SetFilePosition(fptr_specific, "//== SelectData_iv ==//"))
   {
      if ( fscanf(fptr_specific, " %d ", &nvec) != 1)
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== SelectData_iv ==// in the model-dependent input data file datainp_TAG.prn");
      dsgemodel_ps->SelectData_iv = CreateVector_int(nvec);
      if ( !ReadVector_int(fptr_specific, dsgemodel_ps->SelectData_iv) )
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data vector after the first row below the line //== SelectData_iv ==// in the model-dependent input data file datainp_TAG.prn");
      dsgemodel_ps->SelectData_iv->flag = V_DEF;
   }
   else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== SelectData_iv ==// in the model-dependent input data file datainp_TAG.prn does not exist");
   //+
   if (dsgemodel_ps->SelectData_iv->n != dsgemodel_ps->ny)
      fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): SelectData_iv->n must equal ny");
   //---
   dsgemodel_ps->dates_dv = CreateDatesVector_lf(4, dsgemodel_ps->BeginEndSample_iv->v[0], dsgemodel_ps->BeginEndSample_iv->v[1], dsgemodel_ps->BeginEndSample_iv->v[2], dsgemodel_ps->BeginEndSample_iv->v[3]);
   dsgemodel_ps->nSample = dsgemodel_ps->dates_dv->n;
   dsgemodel_ps->fss = dsgemodel_ps->nSample - dsgemodel_ps->nlags;



   //=== Reads doubles
   //--- Model-specific setup inputs
   if ( !fn_SetFilePosition(fptr_specific, "//== div ==//") || fscanf(fptr_specific, " %lf ", &dsgemodel_ps->div) != 1 )
      dsgemodel_ps->div = 1.0+1.0e-09;  //Default value if no input value is supplied.
   if ( !fn_SetFilePosition(fptr_specific, "//== scale4logpriordensity ==//") || fscanf(fptr_specific, " %lf ", &dsgemodel_ps->scale4logpriordensity) != 1 )
      dsgemodel_ps->scale4logpriordensity = -1.0;  //Default: computing the log scale by calling tz_GetScaleForLogpriordensity().



   //=== Reads matricies.
   if (fn_SetFilePosition(fptr_rawdata, "//== Datarawe_dm ==//"))
   {
      if ( fscanf(fptr_rawdata, " %d ", &nrows) != 1)
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== Datarawe_dm ==// in the raw data file dataraw_*.prn");
      if ( fscanf(fptr_rawdata, " %d ", &ncolse) != 1)
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the second integer in the first row below the line //== Datarawe_dm ==// in the raw data file dataraw_*.prn");
      Datarawe_dm = CreateMatrix_lf(nrows, ncolse);
      if ( !ReadMatrix_lf(fptr_rawdata, Datarawe_dm) )
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== Datarawe_dm ==// in the raw data file dataraw_*.prn");
      Datarawe_dm->flag = M_GE;
      dsgemodel_ps->nData = Datarawe_dm->nrows;
      //Incorrect from now on. if ((ncolse-1) != dsgemodel_ps->ny )
      //   fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): Make sure that ny and the columns of Datarawe_dm\n"
      //                   "       is the same for the two files pointed by fptr_rawdata and fptr_common");
      //---
      dsgemodel_ps->Dataraw_dm = CreateMatrix_lf(nrows, ncolse-1);
      dsgemodel_ps->Datarawtran_dm = CreateMatrix_lf(ncolse-1, nrows);
      CopySubmatrix0(dsgemodel_ps->Dataraw_dm, Datarawe_dm, 0, 1, Datarawe_dm->nrows, Datarawe_dm->ncols-1);
      TransposeRegular(dsgemodel_ps->Datarawtran_dm, dsgemodel_ps->Dataraw_dm);
   }
   else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== Datarawe_dm ==// in the raw data file dataraw_*.prn does not exist");
   //--- Getting the sample data for estimation.
   dsgemodel_ps->SampleData_dm = CreateMatrix_lf(dsgemodel_ps->nSample, dsgemodel_ps->SelectData_iv->n);
   dsgemodel_ps->SampleDatatran_dm = CreateMatrix_lf(dsgemodel_ps->SelectData_iv->n, dsgemodel_ps->nSample);
   loc_row_samplebegin = fn_locofyearqm(4, dsgemodel_ps->BeginEndData_iv->v[0], dsgemodel_ps->BeginEndData_iv->v[1], dsgemodel_ps->BeginEndSample_iv->v[0], dsgemodel_ps->BeginEndSample_iv->v[1]);
   CopySubmatrixSelectedCols0(dsgemodel_ps->SampleData_dm, Datarawe_dm, loc_row_samplebegin, dsgemodel_ps->nSample, dsgemodel_ps->SelectData_iv);
   TransposeRegular(dsgemodel_ps->SampleDatatran_dm, dsgemodel_ps->SampleData_dm);
   for (i_=dsgemodel_ps->SampleData_dm->nrows*dsgemodel_ps->SampleData_dm->ncols-1; i_>=0; i_--)
      if (dsgemodel_ps->SampleData_dm->M[i_]<-1.0E+299)
         fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the sample data contains NaNs or -1.0E+301");
   if ((dsgemodel_ps->nnans = dsgemodel_ps->nData - dsgemodel_ps->nSample)<0)
      fn_DisplayError("*_comfuns.c/CreateTSdsgemodel(): making sure BeginEndData contains BeginEndSample specified in the common input data file datainp_common*.prn and model-dependent input data file datainp_TAG.prn");
   if (1)
   {
      fprintf(tzGetDebugFile(), "SampleData_dm:\n");
      WriteMatrix(tzGetDebugFile(), dsgemodel_ps->SampleData_dm, " %10.5f ");
   }


   //=== Reads vectors.
   //--- Model-specific setup inputs
   if (indxStartValuesForMin==1)
   {
      if (fn_SetFilePosition(fptr_specific, "//== xphi_dv ==//"))
      {
         if (fscanf(fptr_specific, " %d ", &nvec) != 1)
            fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== xphi_dv ==// in the input data file datainp_TAG.prn");
         dsgemodel_ps->xphi_dv = CreateVector_lf(nvec);
         if ( !ReadVector_lf(fptr_specific, dsgemodel_ps->xphi_dv) )
            fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== xphi_dv ==// in the input data file datainp_TAG.prn");
         dsgemodel_ps->xphi_dv->flag = V_DEF;
      }
      else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== xphi_dv ==// in the input data file datainp_TAG.prn does not exist");
   }
   else
   {
      if (fn_SetFilePosition(fptr_specific, "//== xphi_dv ==//"))
      {
         if ( fscanf(fptr_specific, " %d ", &nvec) != 1)
            fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== xphi_dv ==// in the input data file datainp_TAG.prn");
         dsgemodel_ps->xphi_dv = CreateVector_lf(nvec);
      }
      else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== xphi_dv ==// in the input data file datainp_TAG.prn does not exist");
   }
   //--- Common setup inputs
   //if (dsgemodel_ps->Use_IMSL_opt)
   //{
   // if (dsgemodel_ps->imsl_ncons)  //If there are non-simple constraints.
   // {
   //    if (fn_SetFilePosition(fptr_common, "//== imsl_lh_vals_dv ==//"))
   //    {
   //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== imsl_lh_vals_dv ==// in the input data file");
   //       dsgemodel_ps->imsl_lh_vals_dv = CreateVector_lf(nvec);
   //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->imsl_lh_vals_dv) )
   //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== imsl_lh_vals_dv ==// in the input data file");
   //       dsgemodel_ps->imsl_lh_vals_dv->flag = V_DEF;
   //    }
   //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== imsl_lh_vals_dv ==// in the input data file does not exist");
   //    //+
   //    if (fn_SetFilePosition(fptr_common, "//== imsl_rh_vals_dv ==//"))
   //    {
   //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== imsl_rh_vals_dv ==// in the input data file");
   //       dsgemodel_ps->imsl_rh_vals_dv = CreateVector_lf(nvec);
   //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->imsl_rh_vals_dv) )
   //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== imsl_rh_vals_dv ==// in the input data file");
   //       dsgemodel_ps->imsl_rh_vals_dv->flag = V_DEF;
   //    }
   //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== imsl_rh_vals_dv ==// in the input data file does not exist");
   // }
   // //--- Model-specific setup inputs
   // // if (dsgemodel_ps->indx_hard_restrict==0) //hard0_ tag: default model (persistence for all exogenous shocks)
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard0_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard0_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard0_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard0_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard0_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard0_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard0_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard0_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard0_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard0_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard0_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard0_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard0_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard0_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==1) //hard1_ tag: no persistence on the parameters of price and wage markup procecess.
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard1_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard1_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard1_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard1_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard1_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard1_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard1_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard1_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard1_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard1_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard1_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard1_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard1_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard1_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==2) //hard2_ tag: no persistence on all exogenous processes.
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard2_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard2_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard2_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard2_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard2_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard2_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard2_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard2_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard2_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard2_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard2_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard2_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard2_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard2_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==3) //hard3_
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard3_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard3_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard3_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard3_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard3_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard3_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard3_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard3_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard3_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard3_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard3_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard3_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard3_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard3_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==103) //hard103_
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard103_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard103_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard103_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard103_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard103_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard103_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard103_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard103_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard103_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard103_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard103_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard103_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard103_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard103_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==203) //hard203_
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard203_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard203_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard203_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard203_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard203_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard203_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard203_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard203_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard203_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard203_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard203_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard203_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard203_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard203_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==4) //hard4_ tag: no persistence on all exogenous processes.
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard4_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard4_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard4_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard4_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard4_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard4_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard4_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard4_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard4_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard4_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard4_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard4_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard4_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard4_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==5) //hard5_ tag
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard5_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard5_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard5_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard5_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard5_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard5_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard5_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard5_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard5_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard5_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard5_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard5_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard5_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard5_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==105) //hard105_ tag
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard105_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard105_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard105_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard105_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard105_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard105_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard105_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard105_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard105_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard105_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard105_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard105_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard105_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard105_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==205) //hard205_ tag
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard205_simple_lowvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard205_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard205_simple_lowvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard205_simple_lowvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard205_simple_lowvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard205_simple_lowvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard205_simple_lowvals_dv ==// in the input data file does not exist");
   // //    //+
   // //    if (fn_SetFilePosition(fptr_common, "//== hard205_simple_highvals_dv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard205_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard205_simple_highvals_dv = CreateVector_lf(nvec);
   // //       if ( !ReadVector_lf(fptr_common, dsgemodel_ps->hard205_simple_highvals_dv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard205_simple_highvals_dv ==// in the input data file");
   // //       dsgemodel_ps->hard205_simple_highvals_dv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard205_simple_highvals_dv ==// in the input data file does not exist");
   // // }
   //}


   //=== Reads integer vectors.
   //--- Common setup inputs
   //if (dsgemodel_ps->Use_IMSL_opt)
   //{
   // if (dsgemodel_ps->imsl_ncons)  //If there are non-simple constraints.
   // {
   //    if (fn_SetFilePosition(fptr_common, "//== imsl_lh_locs_iv ==//"))
   //    {
   //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== imsl_lh_locs_iv ==// in the input data file");
   //       dsgemodel_ps->imsl_lh_locs_iv = CreateVector_int(nvec);
   //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->imsl_lh_locs_iv) )
   //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== imsl_lh_locs_iv ==// in the input data file");
   //       dsgemodel_ps->imsl_lh_locs_iv->flag = V_DEF;
   //    }
   //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== imsl_lh_locs_iv ==// in the input data file does not exist");
   // }
   // //--- Model-specific setup inputs
   // // if (dsgemodel_ps->indx_hard_restrict==0) //hard0_ tag: default model (persistence for all exogenous shocks)
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard0_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard0_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard0_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard0_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard0_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard0_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard0_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==1) //hard1_ tag: no persistence on the parameters of price and wage markup procecess.
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard1_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard1_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard1_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard1_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard1_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard1_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard1_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==2) //hard2_ tag: no persistence on all exogenous processes.
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard2_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard2_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard2_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard2_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard2_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard2_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard2_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==3) //hard3_ tag: no persistence on all exogenous processes.
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard3_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard3_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard3_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard3_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard3_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard3_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard3_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==103) //hard103_ tag:
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard103_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard103_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard103_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard103_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard103_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard103_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard103_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==203) //hard203_ tag:
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard203_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard203_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard203_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard203_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard203_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard203_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard203_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==4) //hard4_ tag: no persistence on all exogenous processes.
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard4_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard4_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard4_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard4_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard4_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard4_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard4_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==5) //hard5_ tag
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard5_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard5_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard5_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard5_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard5_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard5_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard5_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==105) //hard105_ tag
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard105_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard105_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard105_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard105_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard105_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard105_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard105_simple_locs_iv ==// in the input data file does not exist");
   // // }
   // // else if (dsgemodel_ps->indx_hard_restrict==205) //hard205_ tag
   // // {
   // //    if (fn_SetFilePosition(fptr_common, "//== hard205_simple_locs_iv ==//"))
   // //    {
   // //       if ( fscanf(fptr_common, " %d ", &nvec) != 1)
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the first integer in the first row below the line //== hard205_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard205_simple_locs_iv = CreateVector_int(nvec);
   // //       if ( !ReadVector_int(fptr_common, dsgemodel_ps->hard205_simple_locs_iv) )
   // //          fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): check the data matrix or vector after the first row below the line //== hard205_simple_locs_iv ==// in the input data file");
   // //       dsgemodel_ps->hard205_simple_locs_iv->flag = V_DEF;
   // //    }
   // //    else  fn_DisplayError(".../*_comfuns.c/CreateTSdsgemodel(): the line with //== hard205_simple_locs_iv ==// in the input data file does not exist");
   // // }
   //}



   //--------- Allocates memory, creates structures, and assign values. ---------
   dsgemodel_ps->dw_ProbStates_dv = CreateVector_lf(dsgemodel_ps->nSample+1);
   dsgemodel_ps->gensys_ps = CreateTSgensys(gensys_sims, dsgemodel_ps->ngensys, dsgemodel_ps->nshocks, dsgemodel_ps->nexperrs, dsgemodel_ps->div);
   //---
   dsgemodel_ps->priorhyperpars_ps = CreateTSpriorhyperpars();
   //--
   dsgemodel_ps->strshks_dm = CreateMatrix_lf(dsgemodel_ps->ny, dsgemodel_ps->fss);
   //--
   dsgemodel_ps->nadded = dsgemodel_ps->nz - dsgemodel_ps->nzbase;

   //====== Markov-switching Kalman filter ======
   dsgemodel_ps->kalfilmsinputs_1stapp_ps = NULL;


   //===
   DestroyMatrix_lf(Datarawe_dm);  //(fss+lags)-by-(1+ny)
   dsgemodel_ps->Dataraw_dm = DestroyMatrix_lf(dsgemodel_ps->Dataraw_dm);  //Reset to NULL -- Important.
   dsgemodel_ps->Datarawtran_dm = DestroyMatrix_lf(dsgemodel_ps->Datarawtran_dm);  //Reset to NULL -- Important.
   dsgemodel_ps->BeginEndData_iv = DestroyVector_int(dsgemodel_ps->BeginEndData_iv);  //Reset to NULL -- Important.


   return (dsgemodel_ps);
}
//===
void InitializeCreateTSdsgemodelFromStateModel(struct TStateModel_tag *smodel_ps)
{
   int _i;
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   //===
   TSivector *rows_iv = NULL;
   TSivector *cols_iv = NULL;


   dsgemodel_ps->nfreepars_tot = dsgemodel_ps->xphi_dv->n + NumberFreeParametersQ(smodel_ps);
   dsgemodel_ps->Omega_dm = CreateMatrix_lf(dsgemodel_ps->nfreepars_tot, dsgemodel_ps->nfreepars_tot);
   dsgemodel_ps->Hessian_dm = CreateMatrix_lf(dsgemodel_ps->nfreepars_tot, dsgemodel_ps->nfreepars_tot);

   //dsgemodel_ps->kalfilmsinputs_1stapp_ps = CreateTSkalfilmsinputs_1stapp(dsgemodel_ps->ny, dsgemodel_ps->nz, smodel_ps->sv->nstates, dsgemodel_ps->nSample);
   dsgemodel_ps->kalfilmsinputs_1stapp_ps = CreateTSkalfilmsinputs_1stapp2(dsgemodel_ps->ny, dsgemodel_ps->nz, dsgemodel_ps->nu, dsgemodel_ps->ne, smodel_ps->sv->nstates, dsgemodel_ps->nSample);

   dsgemodel_ps->PostProbS_dm = CreateMatrix_lf(dsgemodel_ps->fss+1, smodel_ps->sv->nstates);  //(fss+1)-by-smodel_ps->sv->nstates: marginal posterior probabilities of states

   rows_iv = CreateVector_int(smodel_ps->sv->n_state_variables);
   cols_iv = CreateVector_int(smodel_ps->sv->n_state_variables);
   for (_i=smodel_ps->sv->n_state_variables-1; _i>=0; _i--)
   {
      rows_iv->v[_i] = dsgemodel_ps->fss+1;
      cols_iv->v[_i] = smodel_ps->sv->state_variable[_i]->nbasestates;
   }
   dsgemodel_ps->BaseProbS_dc = CreateCell_lf(rows_iv, cols_iv);
   DestroyVector_int(rows_iv);
   DestroyVector_int(cols_iv);
}
//---
struct TSdsgemodel_tag *DestroyTSdsgemodel(struct TSdsgemodel_tag *dsgemodel_ps)
{
   if (dsgemodel_ps)
   {
      DestroyVector_c(dsgemodel_ps->modeltag_cv);
      //
      DestroyMatrix_lf(dsgemodel_ps->Dataraw_dm);  //(fss+lags)-by-ny
      DestroyMatrix_lf(dsgemodel_ps->Datarawtran_dm);
      DestroyVector_int(dsgemodel_ps->BeginEndData_iv);
      //
      DestroyMatrix_lf(dsgemodel_ps->SampleData_dm);
      DestroyMatrix_lf(dsgemodel_ps->SampleDatatran_dm);
      DestroyVector_int(dsgemodel_ps->BeginEndSample_iv);
      DestroyVector_int(dsgemodel_ps->SelectData_iv);
      DestroyVector_lf(dsgemodel_ps->dates_dv);
      //
      DestroyVector_lf(dsgemodel_ps->xphi_dv);
      DestroyMatrix_lf(dsgemodel_ps->Omega_dm);
      DestroyMatrix_lf(dsgemodel_ps->Hessian_dm);

      //=== The order matters!
      DestroyTSDWRegimeControlVariables(dsgemodel_ps->DWRegimeControlVariables_ps);
      DestroyTSfinalpars(dsgemodel_ps->finalpars_ps);
      DestroyTSgensys(dsgemodel_ps->gensys_ps);
      DestroyTSpriorhyperpars(dsgemodel_ps->priorhyperpars_ps);
      //+ IMSL arguments
      DestroyVector_int(dsgemodel_ps->imsl_lh_locs_iv);
      DestroyVector_lf(dsgemodel_ps->imsl_lh_vals_dv);
      DestroyVector_lf(dsgemodel_ps->imsl_rh_vals_dv);
      //+
      DestroyVector_int(dsgemodel_ps->hard7_imsl_lh_locs_iv);
      DestroyVector_lf(dsgemodel_ps->hard7_imsl_lh_vals_dv);
      DestroyVector_lf(dsgemodel_ps->hard7_imsl_rh_vals_dv);
      //+
      DestroyVector_int(dsgemodel_ps->hard0_simple_locs_iv);
      DestroyVector_lf(dsgemodel_ps->hard0_simple_lowvals_dv);
      DestroyVector_lf(dsgemodel_ps->hard0_simple_highvals_dv);
      //+
      DestroyVector_int(dsgemodel_ps->hard1_simple_locs_iv);
      DestroyVector_lf(dsgemodel_ps->hard1_simple_lowvals_dv);
      DestroyVector_lf(dsgemodel_ps->hard1_simple_highvals_dv);
      //+
      DestroyVector_int(dsgemodel_ps->hard2_simple_locs_iv);
      DestroyVector_lf(dsgemodel_ps->hard2_simple_lowvals_dv);
      DestroyVector_lf(dsgemodel_ps->hard2_simple_highvals_dv);
      //+ Markov-switching Kalman filter and others.
      //DestroyTSkalfilmsinputs_1stapp(dsgemodel_ps->kalfilmsinputs_1stapp_ps);
      DestroyTSkalfilmsinputs_1stapp2(dsgemodel_ps->kalfilmsinputs_1stapp_ps);
      //+ Regime-switching output arguments
      DestroyVector_lf(dsgemodel_ps->dw_ProbStates_dv);
      DestroyMatrix_lf(dsgemodel_ps->PostProbS_dm);
      DestroyCell_lf(dsgemodel_ps->BaseProbS_dc);
      DestroyMatrix_lf(dsgemodel_ps->strshks_dm);

      //--- To be destroyed last!
      tzDestroy(dsgemodel_ps);

      return ((struct TSdsgemodel_tag *)NULL);
   }
   else  return (dsgemodel_ps);
}


struct TSDWRegimeControlVariables_tag *CreateTSDWRegimeControlVariables(FILE *fptr_markov, struct TStateModel_tag *smodel_ps)
{
   int ki;
   int successflag;
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   //--- The following order matters.
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = tzMalloc(1, struct TSDWRegimeControlVariables_tag);
   int *n_st_coef = NULL;
   int *cum_total_coef = NULL;
   int **indx_st_coef = NULL;
   int *n_st_var = NULL;
   int *cum_total_var = NULL;
   int **indx_st_var = NULL;

   //--- Getting regime control table.
   successflag = ReadDSGEControllingVariables(&DWRegimeControlVariables_ps->n_constpars,
                                              &n_st_coef, &cum_total_coef, &indx_st_coef, &n_st_var, &cum_total_var, &indx_st_var,
                                              fptr_markov, smodel_ps->sv);
   if (!successflag)
      fn_DisplayError("Fatal error in CreateTSDWRegimeControlVariables() in .../*_comfuns.c:\n"
                      "    failure in creating regime control table from DW's ReadDSGEControllingVariables()");

   DWRegimeControlVariables_ps->n_st_coef = n_st_coef;
   DWRegimeControlVariables_ps->cum_total_coef = cum_total_coef;
   DWRegimeControlVariables_ps->indx_st_coef = indx_st_coef;
   DWRegimeControlVariables_ps->n_st_var = n_st_var;
   DWRegimeControlVariables_ps->cum_total_var = cum_total_var;
   DWRegimeControlVariables_ps->indx_st_var = indx_st_var;


   dsgemodel_ps->indx_tvmodel = 0;
   dsgemodel_ps->indx_const = 1;
   for (ki=dw_DimA(n_st_coef)-1; ki>=0; ki--)
      if (n_st_coef[ki]>1)
      {
         dsgemodel_ps->indx_tvmodel = 1;
         dsgemodel_ps->indx_const = 0;
         break;
      }
   for (ki=dw_DimA(n_st_var)-1; ki>=0; ki--)
      if (n_st_var[ki]>1)
      {
         dsgemodel_ps->indx_const = 0;
         break;
      }

   return (DWRegimeControlVariables_ps);
}
//---
static struct TSDWRegimeControlVariables_tag *DestroyTSDWRegimeControlVariables(struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps)
{
   if (DWRegimeControlVariables_ps)
   {
      //=== The order matters!
      dw_FreeArray(DWRegimeControlVariables_ps->n_st_coef);
      dw_FreeArray(DWRegimeControlVariables_ps->cum_total_coef);
      dw_FreeArray(DWRegimeControlVariables_ps->indx_st_coef);
      dw_FreeArray(DWRegimeControlVariables_ps->n_st_var);
      dw_FreeArray(DWRegimeControlVariables_ps->cum_total_var);
      dw_FreeArray(DWRegimeControlVariables_ps->indx_st_var);
      //--- To be destroyed last!
      tzDestroy(DWRegimeControlVariables_ps);

      return ((struct TSDWRegimeControlVariables_tag *)NULL);
   }
   else  return (DWRegimeControlVariables_ps);
}















//************************************************************
//** Other key functions.
//************************************************************
//------------------------------------
// Checking other obvious errors on dimensions.
//------------------------------------
void FindFatalErrorsOnDimensions(struct TSdsgemodel_tag *dsgemodel_ps)
{
   int nzbase = dsgemodel_ps->nzbase;
   int indx_tvmodel = dsgemodel_ps->indx_tvmodel;
   TSdvector *xphi_dv = dsgemodel_ps->xphi_dv;
   //---
   struct TSfinalpars_tag *finalpars_ps = dsgemodel_ps->finalpars_ps;
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   int *n_st_var = DWRegimeControlVariables_ps->n_st_var;
   int *cum_total_var = DWRegimeControlVariables_ps->cum_total_var;
   //--- Sims's gensys structure.
   struct TSgensys_tag *gensys_ps = dsgemodel_ps->gensys_ps;
   TSdmatrix *A_dm = gensys_ps->G0_dm;
   int neqs = A_dm->nrows;
   TSdmatrix *Psi_dm = gensys_ps->Psi_dm;
   int nshocks = Psi_dm->ncols;


   if (xphi_dv->n != (cum_total_var[finalpars_ps->n_tvvarpars-1] + n_st_var[finalpars_ps->n_tvvarpars-1]))
         fn_DisplayError(".../*_comfuns.c/FindFatalErrorsOnDimensions(): number of `xphi_dv' read in from datainp_*.prn\n"
                         "     must match the total length of control variables in datainp_markov_*prn");

   if (!indx_tvmodel)  //Case where coefficient parameters are not time varying.
   {
      if (nshocks != finalpars_ps->n_tvvarpars)
         fn_DisplayError(".../*_comfuns.c/FindFatalErrorsOnDimensions(): number of `nshocks' read in from datainp_*.prn\n"
                         "     must match number of 'Time varying variance parameters' in datainp_markov_*prn");
      if (neqs != nzbase)
         fn_DisplayError(".../*_comfuns.c/FindFatalErrorsOnDimensions(): number of `ngensys' read in from datainp_*.prn\n"
                         "     must be nzbase in cases where coefficient parameters are constant");
   }
   else  //Time-varying coefficient parameter case.
   {
      if ((nshocks-n_st_coef[0]) != finalpars_ps->n_tvvarpars)
         fn_DisplayError(".../*_comfuns.c/FindFatalErrorsOnDimensions(): number of `nshocks' read in from datainp_*.prn\n"
                         "     must be 'n_tvvarpars + n_st_coef[0]' where n_st_coef[0] are added dummy regime shocks");
   }

   //--- For all cases.
   if ((dsgemodel_ps->nz - dsgemodel_ps->nadded) != nzbase)
      fn_DisplayError(".../*_comfuns.c/FindFatalErrorsOnDimensions(): Dimension of state and measurement equations\n"
                      "     in kalman filter do not match");
}




//------------------------------------
// Running gensys given NEW parameters.
//------------------------------------
static void RunningGensys_allcases(struct TStateModel_tag *smodel_ps)
{
   // The model written in the canonical form:
   //  A x_t = B x_{t-1} + Gamma_u u_t + Gamma_eta eta_t
   // Output format: x_t = G1*x_{t-1} + impact*u_t.

   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   struct TSfinalpars_tag *finalpars_ps = dsgemodel_ps->finalpars_ps;
   //--- Transition probability for coefficient regimes.
   // TSdmatrix *Qc_dm = finalpars_ps->Qc_dm;
   // TSdvector *qc_ergodic_dv = finalpars_ps->qc_ergodic_dv;
   // int nrows4Qc = Qc_dm->nrows;

   //--- Sims's gensys structure.
   struct TSgensys_tag *gensys_ps = dsgemodel_ps->gensys_ps;
   TSdmatrix *A_dm = gensys_ps->G0_dm;
   TSdmatrix *B_dm = gensys_ps->G1_dm;
   TSdmatrix *Psi_dm = gensys_ps->Psi_dm;
   TSdmatrix *Pi_dm = gensys_ps->Pi_dm;
   int indx_tvmodel = dsgemodel_ps->indx_tvmodel;
   int neqs = A_dm->nrows;


   //------------------------------------------------------------------------------------//
   //---------------- Filling the matrices with new parameters. -------------------------//
   //------------------------------------------------------------------------------------//
   //--- Transition matrix for inflation target.
   // sv_ComputeTransitionMatrix(0, smodel_ps->sv, smodel_ps);  //0 mean time 0 (in the time-varying Q, t means time t)
   //          //This line must be used before copying it to any Q in case theta or q or anything has changed.
   // CopyMatrix0(Qc_dm, smodel_ps->sv->state_variable[0]->baseQ);  //$$##  coefficient regime is ALWAYS specified first.
   // Qc_dm->flag = M_GE;
   // //--- Getting ergodic distribution for coefficient regimes and average inflation.
   // ergodicp(qc_ergodic_dv, Qc_dm);
   // finalpars_ps->log_pistarbar = VectorDotVector(qc_ergodic_dv, log_pistar_dv);


   //--- Constant terms and div.
   InitializeConstantVector_lf(gensys_ps->c0_dv, 0.0);
   gensys_ps->div = dsgemodel_ps->div;

   //=== Input A_dm, B_dm, Psi_dm, and Pi_dim from CInputGensysForm_TAG.cexps exported by Mathematica.
   #include CInputGensysForm_TAG

   //--- Running gensys
   gensys_ps->gensys(gensys_ps, (void *)NULL);
}



//------------------------------------
// For constant and Markov-switching model:
//   (0) must run RunningGensys_allcases() before calling this function, because
//         qc_ergodic_dv is required in this function for the MS model.
//   (1) refreshes the kalman-filter structure with new parameter values;
//   (2) puts shock variances back to the reduced-form solution before forming the LH.
//------------------------------------
static void Refresh_kalfilms_allcases(struct TStateModel_tag *smodel_ps)
{
   //See pp.45-48, LWZ Model II Notes
   //bt_dm, Ft_dc, and Vt_dc MANUALLY change with regime-switching features.

   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSgensys_tag *gensys_ps = dsgemodel_ps->gensys_ps;
   struct TSfinalpars_tag *finalpars_ps = dsgemodel_ps->finalpars_ps;
   struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps = dsgemodel_ps->kalfilmsinputs_1stapp_ps;
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   //int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   //--- Working variables.
   //--- Local variables.
   int si;
   static int indx_once = 0;
   //--- Accessed variables.
   int **regime_index = smodel_ps->sv->index;
   int indx_tvmodel = dsgemodel_ps->indx_tvmodel;
   int nst = kalfilmsinputs_1stapp_ps->nst;
   int ny = kalfilmsinputs_1stapp_ps->ny;
   int nz = kalfilmsinputs_1stapp_ps->nz;
   int nu = kalfilmsinputs_1stapp_ps->nu;
   int ne = kalfilmsinputs_1stapp_ps->ne;
        //In general, ne = nshocks.  But if measurement errors are an AR process, then the i.i.d.
        //  shocks to these errors would be grouped into a part of fundamental shocks.  In this
        //  case, ne > nshocks (which is the dimension for gensys).
   int **indx_st_var = DWRegimeControlVariables_ps->indx_st_var;
   //===
   TSdvector *at_dv = CreateVector_lf(ny);
   TSdmatrix *Wnybynz_dm = CreateMatrix_lf(ny, nz);
   //TSdmatrix *Wnzbynz_dm = CreateMatrix_lf(nz, nz);
   //-
   int nzbase = dsgemodel_ps->nzbase; //the row dimension of the 1st gensys block that does not involve regime-switching coefficients.
   int nshocks = dsgemodel_ps->nshocks; //=nshocks: number of fundamental shocks.
   TSdvector *wne_dv = CreateVector_lf(ne);
   TSdmatrix *Wnebyne_dm = CreateMatrix_lf(ne,ne); //ne-by-ne.
   TSdmatrix *G11_2_dm = CreateMatrix_lf(nzbase, nshocks); //G11^2: first upper block of impact according to p.37a in Liu, Waggoner, and Zha's Great Moderation.
   TSdmatrix *Wimpact_enlarged_dm = CreateMatrix_lf(nz, ne); //Including as a submatrix the first upper block of impact according to p.37a in Liu, Waggoner, and Zha's Great Moderation.
   //TSdmatrix *Wimpact_dm = CreateMatrix_lf(nzbase, nshocks);  //Same dimension of the 1st block gensys impact.
   TSdvector bt_sdv;


   bt_sdv.n = nz;
   bt_sdv.flag = V_DEF;

   if (ne < nshocks)  fn_DisplayError("./Refresh_kalfilms_allcases(): Make sure that ne >= nshocks");


   //--- Filling yt_dm.
   CopyMatrix0(kalfilmsinputs_1stapp_ps->yt_dm, dsgemodel_ps->SampleDatatran_dm);
   //--- Filling at_dm and Ht_dc.
   InitializeConstantMatrix_lf(Wnybynz_dm, 0.0);
   //--- Filling nz-by-nz-by-nst Ft_dc (i.e., G1 for gensys solution z_t = G_1 z_{t-1} + ... .
   InitializeConstantMatrix_lf(kalfilmsinputs_1stapp_ps->Ft_dc->C[0], 0.0);
   CopySubmatrix(kalfilmsinputs_1stapp_ps->Ft_dc->C[0], 0, 0, gensys_ps->Theta_dm, 0, 0, nzbase, nzbase);
   //--- Filling nz-by-ne-by-nst Psiet_dc.
   InitializeConstantMatrix_lf(Wimpact_enlarged_dm, 0.0);

   //=== Importing at_dv->v, Wnybynz_dm->M, Wimpact_enlarged_dm (if ne>nshocks),
   //===   and kalfilmsinputs_1stapp_ps->Ft_dc->C[0]->M from the following file.
   #include CInputMeasureMatrices_TAG


   at_dv->flag = V_DEF;
   for (si=nst-1; si>=0; si--)  CopySubvector2matrix(kalfilmsinputs_1stapp_ps->at_dm, 0, si, at_dv, 0, at_dv->n);
   //+
   for (si=nst-1; si>=0; si--)  CopyMatrix0(kalfilmsinputs_1stapp_ps->Ht_dc->C[si], Wnybynz_dm);
                                //Done with all W*_dim.
   //+
   for (si=nst-1; si>0; si--)  CopyMatrix0(kalfilmsinputs_1stapp_ps->Ft_dc->C[si], kalfilmsinputs_1stapp_ps->Ft_dc->C[0]);
                         //Note that si>0, NOT si>=0 because Ft_dc->C[0] has been used already.
                                          //Done with W*_dm.


   if (indx_once==0)
   {
      //--- Filling Psiut_dc and Rt_dc.
      if (!nu)  //No measurement errors, so no value assignment to kalfilmsinputs_1stapp_ps->Psiut_dc.
      {
         InitializeConstantCell_lf(kalfilmsinputs_1stapp_ps->Rt_dc, 0.0);
      }
      else
      {
         //With measurement errors, kalfilmsinputs_1stapp_ps->Psiut_dc will be passed by
         //  #include "Cinput_MeasureEqns.exps"
         //in the above.
         for (si=nst-1; si>=0; si--)
         {
            fn_DisplayError("./Refresh_kalfilms_allcases(): I have not got time to fill Psiut_dc in CInputMeasureMatrices_TAG");
            MatrixTimesSelf(kalfilmsinputs_1stapp_ps->Rt_dc->C[si], 'U', kalfilmsinputs_1stapp_ps->Psiut_dc->C[si], 'N', 1.0, 0.0);
            SUtoGE(kalfilmsinputs_1stapp_ps->Rt_dc->C[si]); //Making it symmetric and general matrix.
         }
      }

      //--- Filling Gt_dc for initial conditions.
      InitializeConstantCell_lf(kalfilmsinputs_1stapp_ps->Gt_dc, 0.0);

      indx_once = 1;   //After this reset, these commands will not be called again, even when this function is called again and again.
   }

   //--- Filling n_z-by-nst bt_dm.
   InitializeConstantMatrix_lf(kalfilmsinputs_1stapp_ps->bt_dm, 0.0);

   //--- Filling Vt_dc.
   for (si=nst-1; si>=0; si--)
   {
      //--- Putting shock standard deviations back into the solution before forming the LH.
      #include CRefresh_Kalman_allcases_TAG

      InitializeConstantMatrix_lf(Wnebyne_dm, 0.0);
      Wnebyne_dm->flag = M_GE | M_UT | M_LT; //Reset the flag to indicate a diagonal matrix.
      tz_DiagMatrix(Wnebyne_dm, wne_dv);
      CopySubmatrix0(G11_2_dm, gensys_ps->Impact_dm, 0, 0, nzbase, nshocks);
      CopySubmatrix0(Wimpact_enlarged_dm, G11_2_dm, 0,0, G11_2_dm->nrows, G11_2_dm->ncols);
             // Note that the part in the last nadded rows and ne-nshocks columns
             //   is imported from the file Cinput_MeasureEqns*.exps.
      MatrixTimesMatrix(kalfilmsinputs_1stapp_ps->Psiet_dc->C[si], Wimpact_enlarged_dm, Wnebyne_dm, 1.0, 0.0, 'N', 'N');
      //MatrixTimesMatrix(Wimpact_dm, G11_2_dm, Wnebyne_dm, 1.0, 0.0, 'N', 'N');
      //---
      //InitializeConstantMatrix_lf(Wnzbynz_dm, 0.0);
      //CopySubmatrix0(Wnzbynz_dm, Wimpact_dm, 0, 0, Wimpact_dm->nrows, Wimpact_dm->ncols);
      //
      //MatrixTimesSelf(kalfilmsinputs_1stapp_ps->Vt_dc->C[si], 'U', Wnzbynz_dm, 'N', 1.0, 0.0);
      MatrixTimesSelf(kalfilmsinputs_1stapp_ps->Vt_dc->C[si], 'U', kalfilmsinputs_1stapp_ps->Psiet_dc->C[si], 'N', 1.0, 0.0);
      SUtoGE(kalfilmsinputs_1stapp_ps->Vt_dc->C[si]); //Making it symmetric and general matrix.
   }
                                          //Done with W*_dm.
   //===
   DestroyVector_lf(at_dv);
   DestroyMatrix_lf(Wnybynz_dm);
   //DestroyMatrix_lf(Wnzbynz_dm);
   //---
   DestroyVector_lf(wne_dv);
   DestroyMatrix_lf(Wnebyne_dm);
   DestroyMatrix_lf(G11_2_dm);
   DestroyMatrix_lf(Wimpact_enlarged_dm);
   //DestroyMatrix_lf(Wimpact_dm);
}



//------------------------------------
// Refreshes everything that needs be computed, given new parameters.
// Reset the refreshed notification variable (dsgemodel_ps->refereshed).
//------------------------------------
int RefreshEverything(struct TStateModel_tag *smodel_ps)
{
   //Return 1: means that we have refreshed everything successfully.
   //       0: fails to refresh everything.

   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps = dsgemodel_ps->kalfilmsinputs_1stapp_ps;

   dsgemodel_ps->refreshed = 1;
   dsgemodel_ps->ValidParameters = 0;  //1: parameters (bar transition matrix elements) are valid (satisfying a priori restrictions),
                                      //0: not valid.
                                      //This has a meaning ONLY when refreshed=1.
   dsgemodel_ps->ValidGensysSolution = 0;  //1: gensys has a (MSV or unique) solution,
                                          //0: no such a solution,
                                          //This has a meaning ONLY when refreshed=1.
   dsgemodel_ps->ValidInitial_z10_P10 = 0; //1: InitializeKalman_z10_P10() gives valid initial values z_{1|0} and P_{1|0}.
                                          //0: no such valid initials.
                                         //This has a meaning ONLY when refreshed=1.
   kalfilmsinputs_1stapp_ps->ztm1_track = -1;
            //t-1 = -1:      no initial conditions z_{1|0} and P_{1|0} has been computed yet, but will be using InitializeKalman_z10_P10(),
            //t-1 >= 0:T-1:  z_{t|t-1} and P_{t|t-1} are updated up to t-1.
   kalfilmsinputs_1stapp_ps->dtm1_track = -1;
            //t-1 = -1:      no etdata_dc->C[0] or Dtdata_d4->F[0] has been computed yet.
            //t-1 >= 0:T-1:  etdata_dc->C[t-1] and Dtdata_d4->F[t-1] are updated up to t-1.



   //------------------------
   //--- The following functions must be called only ONCE for the integrated-out-over-regimes LH from tbase0=0 to tbase0=T-1.
   //------------------------
   Convertphifreepars2finalpars(dsgemodel_ps);  //Refresh the model's final parameters (bar transition matrix elements)

   //--- Early exit with zero likelihood.
   if (ViolateAPrioriRestrictions(dsgemodel_ps))
   {
      fprintf(tzGetDebugFile(), "\n----- ViolateAPrioriRestrictions() in *_comfuns.c: Restricions violated. -----\n");
      return (0);
   }
   else  dsgemodel_ps->ValidParameters = 1;

   RunningGensys_allcases(smodel_ps);
   //--- Early exit if the solution does not exist.
   if (!dsgemodel_ps->gensys_ps->eu_iv->v[0])
   {
      fprintf(tzGetDebugFile(), "\n----- No solution (no existence) from gensys. -----\n");
      return (0);
   }
   else  dsgemodel_ps->ValidGensysSolution = 1;

   //--- Refresh_kalfilms_allcases() must be called first before InitializeKalman_z10_P10() is called.
   Refresh_kalfilms_allcases(smodel_ps);
   if (!InitializeKalman_z10_P10(kalfilmsinputs_1stapp_ps, (TSdmatrix *)NULL, (TSdcell *)NULL))  return (0);
   else  dsgemodel_ps->ValidInitial_z10_P10 = 1;

   return (1);  //Reset to 1 so that we don't compute these needlessly.
}



//------------------------------------
// log conditional likelihood at time t for both constant and Markov-switching models.
// To be used by DW's Markov-switching code.
//------------------------------------
double logTimetCondLH(int s_t, int inpt, struct TStateModel_tag *smodel_ps)
{
   //s_t: base-0: deals with the cross-section values of s_t at the time t.
   //inpt: base-1 in the sense that inpt>=1 to deal with the time series situation where S_T is (T+1)-by-1 and Y_T is T+nlags_max-by-1.
   //      The 1st element for S_T is S_T[1] while S_T[0] is s_0.  The same for (T+1)-by-1 gbeta_dv and nlcoefs-by-(T+1) galpha_dm.
   //      The 1st element for Y_T, however, is Y_T[nlags_max+1-1].
   //See (42.3) on p.42 in the SWZII NOTES.

   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps = dsgemodel_ps->kalfilmsinputs_1stapp_ps;

   //--- The following calls are absolutely NECESSARY because it means that RefreshEverthing()
   //---   will be called only ONCE for the integrated-out-over-regimes LH from tbase0=0 to tbase0=T-1.
   //---   until the reset functions tz_thetaChanged() and tz_TransitionMatrixChanged() are called.
   //---   DW's code will call tz_thetaChanged() and tz_TransitionMatrixChanged() only when
   //---   parameters are changed (e.g., when inpt=T (reaches the end of the sample)
   //---   and new parameters are given and the forward recursion starts from the beginning again.).
   if (!dsgemodel_ps->refreshed)
   {
      if (!RefreshEverything(smodel_ps))  return (-NEARINFINITY); //Early exit to achieve efficiency.
   }
   else //Early exit to achieve efficiency.
   {
      if (!dsgemodel_ps->ValidParameters || !dsgemodel_ps->ValidGensysSolution || !dsgemodel_ps->ValidInitial_z10_P10)  return (-NEARINFINITY);
   }

   return ( logTimetCondLH_kalfilms_1stapp(s_t, inpt, kalfilmsinputs_1stapp_ps, smodel_ps) );
}




//---
static struct TSpriorhyperpars_tag *DestroyTSpriorhyperpars(struct TSpriorhyperpars_tag *priorhyperpars_ps)
{
   if (priorhyperpars_ps)
   {
      //+  To be destroyed last!
      tzDestroy(priorhyperpars_ps);

      return ((struct TSpriorhyperpars_tag *)NULL);
   }
   else  return (priorhyperpars_ps);
}



//------------------------------------
// One-time call: getting the scale of log prior density for possible hard restrictions on the parameter space.
//------------------------------------
double tz_GetScaleForLogpriordensity(struct TStateModel_tag *smodel_ps)
{
   int drawi;
   double inverse_scale2priordensity;
   //--- Taken out of Waggoner's model.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   TSdvector *xphi_dv = dsgemodel_ps->xphi_dv;


   xphi_dv->flag = V_DEF; //Because it will have values below.
   dsgemodel_ps->nViolates = 0;
   for (drawi=dsgemodel_ps->nDrawsFromPrior-1; drawi>=0; drawi--)
   {
      Draw_xphi_dvFromPrior_Specific(smodel_ps);
      if (!RefreshEverything(smodel_ps))  dsgemodel_ps->nViolates++;
   }


   if ((inverse_scale2priordensity=1.0 - ((double)dsgemodel_ps->nViolates)/((double)dsgemodel_ps->nDrawsFromPrior))<MACHINEZERO)
      return (-NEARINFINITY);
   else
   {
      return (-log(inverse_scale2priordensity));
   }
}


//------------------------------------
// Overall posterior kernal for the constant-parameter model only.
//------------------------------------
double logOverallPosteriorKernal_const(struct TStateModel_tag *smodel_ps, TSdvector *xnew_dv)
{
   //This is only for the constant-parameter case.

   double logvalue;
   //
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;


   //=====================================
   //=== This is a must step to refresh the model-specific parameters at a new point.
   //===                        refresh everything.
   //=====================================
   CopyVector0(dsgemodel_ps->xphi_dv, xnew_dv);
   if (!RefreshEverything(smodel_ps))  return (-NEARINFINITY); //Early exit to achieve efficiency.


   //=====================================
   //=== Computes the overall likelihood and the prior.
   //=====================================
   //--- Overall likelihood pdf.
   logvalue = loglh_const(smodel_ps);   //Constant-parameter case.
   dsgemodel_ps->LHvalatpostpeak = logvalue;  //LH value at the posterior peak.


   //--- OverallPosteriorkernal
   dsgemodel_ps->peaklogpost = ( logvalue += (dsgemodel_ps->logprior = logpriordensity(smodel_ps)) );

   return (logvalue);
}
//------------------------------------
// log likelihood for the constant-parameter model.
// WARNING: making sure that
//      Convertphifreepars2finalpars(dsgemodel_ps);
//   is called before this function is called.
//------------------------------------
static double loglh_const(struct TStateModel_tag *smodel_ps)
{
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   //--- Accessible variables.
   int nSample = dsgemodel_ps->nSample;  //fss+lags.
   int nlags = dsgemodel_ps->nlags;
   //--- Local variables.
   int ti, tp1i;
   double loglh;

   //--- The likelihood using the data from the period nlags+1 to nSample.
   //--- Initial conditions z_{1|0} and P_{1|0} are obtained using lagged data.
   loglh = 0.0;
   for (ti=nlags; ti<nSample; ti++)
   {
      tp1i = ti+1;  //Base-1 timing.
      loglh += logTimetCondLH(0, tp1i, smodel_ps);
   }

   return (loglh);
}




//------------------------------------------------------------
// Functions required by the ThetaRoutines structure or CreateParameters() in switch.c.
//------------------------------------------------------------
int NumberOfFreeModelSpecificParameters(struct TStateModel_tag *smodel_ps)
{
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   return (dsgemodel_ps->xphi_dv->n);
}
//---
void ConvertFreeParameters2ModelSpecificParameters(struct TStateModel_tag *smodel_ps, PRECISION *xopt_pd)
{
   //xopt_pd: the (temporary) vector of free model parameters for optimization or minimization.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   memcpy(dsgemodel_ps->xphi_dv->v, xopt_pd, dsgemodel_ps->xphi_dv->n*sizeof(double));
   Convertphifreepars2finalpars(dsgemodel_ps);
}
//---
void CopyMyFreeParameters2OptimizationParameters(struct TStateModel_tag *smodel_ps, PRECISION *xopt_pd)
{
   //This function allows DW code to have access to dsgemodel_ps->xgphi_dv->v, because in his
   //  generic code, the structure is declared as void *, thus DW cannot have access to these
   //  parameters directly.
   //xopt_pd: the (temporary) vector of free model parameters for optimization or minimization.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   memcpy(xopt_pd, dsgemodel_ps->xphi_dv->v, dsgemodel_ps->xphi_dv->n*sizeof(double));
}
//===
void tz_thetaChanged(struct TStateModel_tag *smodel_ps)
{
   //theta: model parameters bar transition matrix elements.

   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;

   dsgemodel_ps->refreshed = 0;  //Meaning: we need to refresh everying because the parameters haved changed.
}
void tz_TransitionMatrixChanged(struct TStateModel_tag *smodel_ps)
{
   //Free transition matrix parameters.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;

   dsgemodel_ps->refreshed = 0;  //Meaning: we need to refresh everying because the parameters haved changed.
}


//=== New Constrained Optimization Structure related
struct TSModelConstraints_tag *CreateTSModelConstraints(int n, int num_lin, int num_nlin, TFNonLinConstraints *non_lin_con_func)
{
  struct TSModelConstraints_tag *package_modelconstraints_ps = tzMalloc(1, struct TSModelConstraints_tag);
  
  package_modelconstraints_ps->n = n;
  package_modelconstraints_ps->num_lin = num_lin;
  package_modelconstraints_ps->num_nlin = num_nlin;
  
  package_modelconstraints_ps->lb = CreateConstantVector_lf(n, -NPSOLINF);
  package_modelconstraints_ps->ub = CreateConstantVector_lf(n, NPSOLINF);
  
  if (num_lin > 0)  {
    package_modelconstraints_ps->lin_lb = CreateConstantVector_lf(num_lin, -NPSOLINF);
    package_modelconstraints_ps->lin_ub = CreateConstantVector_lf(num_lin, NPSOLINF);
    package_modelconstraints_ps->A = CreateConstantMatrix_lf(num_lin, n, 0);
  }
  else {
    package_modelconstraints_ps->lin_lb = NULL;
    package_modelconstraints_ps->lin_ub = NULL;
    package_modelconstraints_ps->A = NULL;
  }
  
  if (num_nlin<=0)
  {
    package_modelconstraints_ps->nlin_lb = NULL;
    package_modelconstraints_ps->nlin_ub = NULL;
  }
  else
  {
    package_modelconstraints_ps->nlin_lb = CreateConstantVector_lf(num_nlin, -NPSOLINF);
    package_modelconstraints_ps->nlin_ub = CreateConstantVector_lf(num_nlin, NPSOLINF);
  }
  
  package_modelconstraints_ps->non_lin_const_func = non_lin_con_func;
  
  return package_modelconstraints_ps;
}

