/*
 * Check ??? to finish the remaining job.
 * Check $$## for hard-coded lines.
 * Check <<>> for Dan's new code prototypes.
 *
 * Credit friction model.
 * See Zha Notes pp.  and Generate_gensysform_credit4housing.nb.
 *
 * Written by T. Zha, March 2009.
 * Copywright (c) by 2009.
*/

/**
//=== For debugging purpose.
if (1)
{
   double t_loglht;
   t_loglht = -(0.5*ny)*LOG2PI - 0.5*logdeterminant(Dtdata_dm) - 0.5*VectorDotVector(wny_dv, etdata_dv);
   fprintf(tzGetDebugFile(), " %10.5f\n", t_loglht);

   fprintf(tzGetDebugFile(), "%%st=%d, inpt=%d, and sti=%d\n", st, inpt, sti);

   fprintf(tzGetDebugFile(), "wP0_dv:\n");
   WriteVector(tzGetDebugFile(), wP0_dv, " %10.5f ");
   fprintf(tzGetDebugFile(), "Vt_dc->C[sti_v=%d]:\n", sti_v);
   WriteMatrix(tzGetDebugFile(), Vt_dc->C[sti_v], " %10.5f ");

   fflush(tzGetDebugFile());
}
**/

#include <cstring>
#include <sstream>
#include <mpi.h>
#include "dw_rand.h"
#include "CStorageHead.h"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CEquiEnergy_TState.h"
#include "CMetropolis.h"
#include "TaskScheduling.h"

extern "C" 
{
	#include "dsgelinv1_comfuns.h"
	#include "optpackage.h"
	// header files included by HW for equi-energy sampling: start
	#include "dw_parse_cmd.h"
	// header files included by HW for equi-energy sampling: end
}

using namespace std; 

static char MODELTAG[TAGLEN];

//------ For the general minimization problem. -------
static char filename_sp_vec_minproj[TAGLEN];   //This is a string because the memory is on the stack using [256].
       //If csminwel, file name for starting point of vectorized (vec) parameters from the minimization problems with the /c 0 option. Must be allocated big enouch a string area.  Increase the array length if necessary.
       //2 vectors.  1st row: gradient; 2nd row: vectorized parameters.
       //If IMSLconlin, writing is different.

//------ For the csminwel minimization problem. -------
//--- Step 1. Creats a number of csminwel structures for both Markov-switching and constant-parameter models.
static double minobj(struct TSminpack_tag *minpack_ps);  //This function is for the constant-parameter model only.
//--- Step 2.
static void InitializeForMinproblem(int n_args_cl, char **args_cl, struct TSminpack_tag *minpack_ps, char *filename_inpdatasp);
//--- Step 3.
//For the constant-parameter model, run minfinder(minpack_ps);  //Constant-parameter case.
//For the regime-switching model, run minfinder_blockcsminwel(minpack_ps);  //Time-varying case.


//------ For IMSL linearly-constrained minimization. -------
//--- Step 1. Calling CreateTSpackagae_imslconlin() and finishing other setups.
//--- Step 2. Calling InitializeForMinproblem();
//--- Step 3. Calling minfinder_noblockimslconlin();



//--- Other static function declarations.
//static void ftd_const_CounterFactual(FILE *fptr_matlab3, struct TStateModel_tag *smodel_ps);
static void ftd_FprintMatlab2(FILE *fptr_matlab2, struct TStateModel_tag *smodel_ps); //For debugging purpose -- printing out the likelihood at time t and regime st.
static void ftd_FprintMatlab1(FILE *fptr_matlab1, struct TStateModel_tag *smodel_ps);
static void ftd_FprintMatlab(FILE *fptr_matlab, struct TStateModel_tag *smodel_ps); //Probabilities of regimes.
static void ftd_FprintEstResults(FILE *fptr_kernel7x_dv_localmodes, TSminpack *minpack_ps);
static void ftd_FprintOutput(FILE *fptr_output, TSminpack *minpack_ps);


// Equi-energy sampler, master and slave processes: start
// void master_single_thread(int argc, char **argv, CEESParameter &parameter, CStorageHead &storage, int highest_level, int lowest_level, int deposit_frequency, int simulation_length, int burn_in_length, int max_energy_tuning_time);
// parameter cannot be const, because its h0 may be updated
// storage cannot be const, because it is consolidated or cleared at the master node
// tuning at different levels are carried out in parallel

// void master_serial_tuning(int argc, char** argv, CEESParameter &parameter, CStorageHead &storage, int highest_level, int lowest_level, int deposit_frequency, int simulation_length, int burn_in_length, int max_energy_tuning_time);
// during the first time of energy tuning, tuning at different levels
// are carried out in serial
// with level of the lowest temperature starting first, 
// followed by levels of higher temperatures, because
// .theta.variance.i will be used as the starting point of 
// .theta.variance.(i+1) during the first time of energy 
// tuning. In the later tuning iterations, tuning are carried out in
// parallel, because .theta.vairance.i of the previous iteration 
// will be used as the starting point of .theta.variance.i of the 
// current iteration

// void slave_single_thread(int argc, char **argv, CEESParameter &parameter, CStateModel &target, CStorageHead &storage, const gsl_rng *r);
// parameter cannot be const, because every time a slave node receives a message from the master
// it will extract h0 from the received message and then update parameter.h and parameter.t
// storage cannot be const, because its bins are constantly altered for deposition and drawing
// target cannot be const, because its energy-level, energy and temperature bounds will be changed
// physically, and its model will be changed logically.

// Equi-energy sampler, master and slave processes: end

//----------------------------------------------------
//-- Main program begins here.
//----------------------------------------------------
int main(int n_args_cl, char **args_cl)
{
   char *cl_modeltag = NULL,   //Point to the command line model tag.
        *cl_filename = NULL,      //Pointer to the command-line raw data file name or filenamebuffer,
                                  //  which moves around (i.e., NOT fixed unless its memory is allocated and string copied.
        filenamebuffer[TAGLEN];  //Used for multiple file names, so must be allocated big enouch a string area.  Increase the array length if necessary.
        TScvector *filename_inpdatasp_cv = CreateVector_c(TAGLEN);
        FILE *fptr_rawdata = NULL,      //File pointer to raw data input file.
        *fptr_common = NULL,    //File pointer to common setup input file.
        *fptr_specific = NULL,    //File pointer to model-specific setup input file.
        *fptr_kernel7x_dv_localmodes = (FILE *)NULL,  //Writes log posterior kernel and the mode estimates x_dv for each starting point to a row of the file outdataest_localmodes_TAG.prn.
        *fptr_output = NULL,     //Writes final output arguments to a file to be examined by us.
        *fptr_matlab = NULL,   //File pointer to output data file for forecasts to be used by Matlab.
        *fptr_matlab1 = NULL,   //File pointer to output data file for beliefs and forecasts to be used by Matlab.
        *fptr_matlab2 = NULL,   //File pointer to output data file for forecasts to be used by Matlab.
        *fptr_matlab3 = NULL,   //File pointer to output data file that report artificial data.
        *fptr_markov = NULL;   //File pointer to setup input data file for Markov switching models, using Waggoner's code.
   //--- Local variables.
   int drawi;
   int indxFindMLE = 0;   //1: find MLE without a prior, 0: find posterior (with a prior).
   int indxStartValuesForMin;  
         //  0: read in from the last estimated results contained in filename_sp if it exists.
         //  1: starts from the fixed values for xphi_dv, read in from datainp_*.prn.
         //  2: randomly or arbitarily selects the initial starting values for the MLE or posterior estimate.
   int indxInitializeTransitionMatrix;
   //--- My model structure.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)NULL;
   //--- Waggoner's Markov switching package.
   struct TMarkovStateVariable_tag *sv_ps = (struct TMarkovStateVariable_tag *)NULL;
   ThetaRoutines *sroutines_ps = (ThetaRoutines *)NULL;
   struct TStateModel_tag *smodel_ps = (struct TStateModel_tag *)NULL;
   struct TStateModel_tag *dw_model=(TStateModel*)NULL;
   //--- General (csminwel) minimization for constant-parameter.
   struct TSetc_minproj_tag *etc_minproj_ps = NULL;
   struct TSminpack_tag *minpack_ps = NULL;
   //--- Blockwise (csminwel) minimization for regime-switching model.
   struct TSargs_blockcsminwel_tag *args_blockcsminwel_ps = NULL;
   #if defined(_IMSLOPTPACKAGE_)
   //--- IMSL minimization.
   int ki;
   struct TSpackage_imslconlin_tag *package_imslconlin_ps = NULL;
   #endif
   //--- Regime-switching output arguments.
   int si, svi, si_subgrand, s_t_i, basei;
   TSdvector postprobs_sdv, postprobs2_sdv;
	
	
   //-----------------
   // Reads from the command line the user-specified input file and the most-often-used integer arguments such as sample size.
   //-----------------
   cl_modeltag = fn_ParseCommandLine_String(n_args_cl,args_cl,'t',(char *)NULL);  // Tag for different models.
   if (!cl_modeltag)  //Must print this out on the screen because tzGetDebugFile() is not yet defined, which needs cl_modeltag to be defined anyway.
   {
      printf("\n &&&&&&&&&&&&&&&&&&&&&& \n");
      printf("Fatal error: .../main(): No model tag cl_modeltag is specified in the command line!");
      fflush(stdout);
   }
   else  strcpy(MODELTAG, cl_modeltag); //MODELTAG is used in ftd_FprintOutput().
   indxStartValuesForMin = fn_ParseCommandLine_Integer(n_args_cl,args_cl,'c',1);
   sprintf(filename_sp_vec_minproj, "outdatasp_min_%s.prn", cl_modeltag);
   //+
   sprintf(filenamebuffer, "dataraw.prn");
   cl_filename = fn_ParseCommandLine_String(n_args_cl,args_cl,'r',filenamebuffer);  //Raw data input file.
   fptr_rawdata = tzFopen(cl_filename,"r");
   //+
   sprintf(filenamebuffer, "datainp_common.prn");
   cl_filename = fn_ParseCommandLine_String(n_args_cl,args_cl,'i',filenamebuffer);  //Common setup input data file.
   fptr_common = tzFopen(cl_filename,"r");
   //+
   sprintf(filenamebuffer, "datainp_%s.prn", cl_modeltag);
   cl_filename = fn_ParseCommandLine_String(n_args_cl,args_cl,'s',filenamebuffer);  //Model-specific setupt input data file.
   fptr_specific = tzFopen(cl_filename,"r");
   //+
   sprintf(filenamebuffer, "datainp_markov_%s.prn", cl_modeltag);
   cl_filename = fn_ParseCommandLine_String(n_args_cl,args_cl,'m',filenamebuffer);  //Markov-switching setup input data file.
   fptr_markov = tzFopen(cl_filename,"r");
   //+
   sprintf(filename_inpdatasp_cv->v, "inpdatasp_min_%s.prn", cl_modeltag);
   filename_inpdatasp_cv->flag = V_DEF;
   //--- Output data files.
   //--- Output data files.
   sprintf(filenamebuffer, "outdata_debug_%s.prn", cl_modeltag);
   tzOpenDebugFile(filenamebuffer);  //Debug output file.
   //+
   sprintf(filenamebuffer, "outdata_optinterm_%s.prn", cl_modeltag);
   tzOpenOptIntermFile(filenamebuffer);  //Debug output file.
   //+
   sprintf(filenamebuffer, "outdataout_%s.prn", cl_modeltag);
   fptr_output = tzFopen(filenamebuffer,"w");  //Final output file.
   //+
   sprintf(filenamebuffer, "outdataest_localmodes_%s.prn", cl_modeltag);
   fptr_kernel7x_dv_localmodes = tzFopen(filenamebuffer,"w");  //Output file for log posterior kernel and estimates at a local mode.
   //+
   sprintf(filenamebuffer, "outdatainp_matlab_%s.prn", cl_modeltag);
   fptr_matlab = tzFopen(filenamebuffer, "w");
   //+
   sprintf(filenamebuffer, "outdatainp_matlab1_%s.prn", cl_modeltag);
   fptr_matlab1 = tzFopen(filenamebuffer, "w");
   //+
   sprintf(filenamebuffer, "outdatainp_matlab2_%s.prn", cl_modeltag);
   fptr_matlab2 = tzFopen(filenamebuffer, "w");
   //+
   sprintf(filenamebuffer, "outdatainp_matlab3_%s.prn", cl_modeltag);
   fptr_matlab3 = tzFopen(filenamebuffer, "w");


   //----------------------------------------------
   //--- Memory allocation and structure creation.
   //--- The order matters!
   //----------------------------------------------
   //====== Creating my model structure (including reading the data and parameter vlaues). ======
   //====== Checking the model type.  ======
   dsgemodel_ps = CreateTSdsgemodel(fptr_rawdata, fptr_common, fptr_specific, fptr_markov, indxStartValuesForMin);
   if (dsgemodel_ps->modeltag_cv->n > TAGLEN)
   {
      printf("\n &&&&&&&&&&&&&&&&&&&&&& \n");
      printf("Fatal error: .../main(): The length of cl_modeltag must be less than dsgemodel_ps->modeltag_cv->n!");
      fflush(stdout);
   }
   strncpy(dsgemodel_ps->modeltag_cv->v, cl_modeltag, dsgemodel_ps->modeltag_cv->n-2);  //-2 to be save.  I think -1 would be enough to allow the terminating null character.
   dsgemodel_ps->modeltag_cv->flag = V_DEF;
   //+
   // if (dsgemodel_ps->indx_hard_restrict==0) //0: default model (persistence for all exogenous shocks)
   // {
   //    if (strncmp("mod0", cl_modeltag,  4))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 0 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==1) //1: no persistence on the parameters of price and wage markup procecess.
   // {
   //    if (strncmp("mod1", cl_modeltag,  4))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 1 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==2) //2: no persistence on all exogenous processes.
   // {
   //    if (strncmp("mod2", cl_modeltag,  4))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 2 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==3)
   // {
   //    if (strncmp("mod3", cl_modeltag,  4))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 2 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==103)
   // {
   //    if (strncmp("modb103", cl_modeltag,  7))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 2 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==203)
   // {
   //    if (strncmp("modc203", cl_modeltag,  7))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 2 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==4)
   // {
   //    if (strncmp("mod4", cl_modeltag,  4))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 2 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==5)
   // {
   //    if (strncmp("mod5", cl_modeltag,  4))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 2 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==105)
   // {
   //    if (strncmp("modb105", cl_modeltag,  7))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 2 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else if (dsgemodel_ps->indx_hard_restrict==205)
   // {
   //    if (strncmp("modc205", cl_modeltag,  7))  fn_DisplayError("main(): indx_hard_restrict in datainp_common.prn must be 2 or\n"
   //                                                           "\t its value must match mod# in the datainp_mod#_*.prn");
   // }
   // else  fn_DisplayError("main():  indx_hard_restrict in datainp_common.prn is NOT properly specified");


   //====== Waggoner's Markov switching variables. ======
   sv_ps = CreateMarkovStateVariable_File(fptr_markov, (char *)NULL);
            //In this case, fptr_markov points to datainp_markov_const.prn, which can be found in D:\ZhaData\CommonFiles\C_Examples_DebugTips\DW_MarkovInputFiles.
   sroutines_ps = CreateThetaRoutines_empty();  //_empty means all the functions are set to NULL and then are supplied below.
   sroutines_ps->pLogConditionalLikelihood = logTimetCondLH;                                     //User's: logTimetCondLH
   sroutines_ps->pLogPrior = logpriordensity;                                            //User's: pLogPrior
   sroutines_ps->pNumberFreeParametersTheta = NumberOfFreeModelSpecificParameters;               //User's: NumberOfFreeModelSpecificParameters,
   sroutines_ps->pConvertFreeParametersToTheta = ConvertFreeParameters2ModelSpecificParameters;  //User's: ConvertFreeParameters2ModelSpecificParameters,
   sroutines_ps->pConvertThetaToFreeParameters = CopyMyFreeParameters2OptimizationParameters;  //User's: CopyMyFreeParameters2OptimizationParameters,
   sroutines_ps->pThetaChanged = tz_thetaChanged;                                                   //User's: Notification routine (need to refresh everything given new parameters?)
   sroutines_ps->pTransitionMatrixParametersChanged = tz_TransitionMatrixChanged;   //User's: Notification routine (need to refresh everything given new parameters?)
   smodel_ps = CreateStateModel(dsgemodel_ps->nSample, sv_ps, sroutines_ps, (void *)dsgemodel_ps);
   smodel_ps->fobs = dsgemodel_ps->nlags + 1; //meaning: we compute the LH from nlags+1 to ndata (i.e., T).
   //--- Optional.
   if (1)
   {
      if (fn_SetFilePosition(fptr_markov, "//== indxInitializeTransitionMatrix ==//"))
         if ((fscanf(fptr_markov, " %d ", &indxInitializeTransitionMatrix) == 1) && indxInitializeTransitionMatrix)
            if(!ReadBaseTransitionMatrices(fptr_markov, (char*)NULL, "Initial: ", smodel_ps))  //Waggoner's function.
               fn_DisplayError(".../main.c: Supply initial transition matrices in datainp_markov*.prn");
   }
   //--- Finishing up the creation of TSdsgemodel -- the order matters!
   InitializeCreateTSdsgemodelFromStateModel(smodel_ps);
   dsgemodel_ps->DWRegimeControlVariables_ps = CreateTSDWRegimeControlVariables(fptr_markov, smodel_ps);
   dsgemodel_ps->finalpars_ps = CreateTSfinalpars(fptr_common, dsgemodel_ps);
   FindFatalErrorsOnDimensions(dsgemodel_ps); //Essential step!


   //====== csminwel minimization setup (step 1). ======
   args_blockcsminwel_ps = CreateTSargs_blockcsminwel(fptr_common);
                //Blockwise (csminwel) minimization arguments, reading convergence criteria or using default values if fptr_common is set to NULL.
                //fptr_common contains parameters for both constant-parameter and Markov-switching models.
   etc_minproj_ps = CreateTSetc_minproj(&smodel_ps, (TFDestroyTStateModel *)NULL, &args_blockcsminwel_ps, DestroyTSargs_blockcsminwel);
                    //Taking convergence criteria and my model structure smodel_ps into minpack_ps.
   minpack_ps = CreateTSminpack((TFminobj *)minobj, (void **)&etc_minproj_ps, (TFmindestroy_etcproject *)NULL, (TFmingrad *)NULL,
                                filename_sp_vec_minproj,
                                dsgemodel_ps->nfreepars_tot,
                                MIN_CSMINWEL);
             //minobj is for the constant-parameter model only in which case, NumberFreeParametersQ(smodel_ps) will be 0.

   //========= IMSL minimization setup (step 1). =========
   //=== $$##: may need to change indx_hard_restrict for the following if loop.
   #if defined(_IMSLOPTPACKAGE_)
   if (dsgemodel_ps->Use_IMSL_opt)
   {
      package_imslconlin_ps = CreateTSpackagae_imslconlin(dsgemodel_ps->nfreepars_tot, dsgemodel_ps->imsl_neqs, dsgemodel_ps->imsl_ncons);
      if (dsgemodel_ps->imsl_ncons)  //If there are non-simple linear constraints.
      {
         if (dsgemodel_ps->imsl_lh_locs_iv->n != dsgemodel_ps->imsl_lh_vals_dv->n)
            fn_DisplayError("main.c: dimensions of imsl_lh_locs_iv and imsl_lh_vals_dv must match");
         for (ki=dsgemodel_ps->imsl_lh_locs_iv->n-1; ki>=0; ki--)
            package_imslconlin_ps->lh_coefs_dv->v[dsgemodel_ps->imsl_lh_locs_iv->v[ki]] = dsgemodel_ps->imsl_lh_vals_dv->v[ki];
         package_imslconlin_ps->lh_coefs_dv->flag = V_DEF;
         //--- Filling in linear constraints rh_constraints_dv.
         if (dsgemodel_ps->imsl_rh_vals_dv->n != dsgemodel_ps->imsl_ncons)
            fn_DisplayError("main.c: # of constraints must match # of rhs values in linear constraints");
         CopyVector0(package_imslconlin_ps->rh_constraints_dv, dsgemodel_ps->imsl_rh_vals_dv);
      }

      //========== Simple bounds. =================
      //--- For free parameters in transition matrices (the last few elments).
      //--- AUTOMATED! i.e., the program will NOT enter this loop for the constant-parameter model.
      //for (ki=dsgemodel_ps->xphi_dv->n; ki<dsgemodel_ps->nfreepars_tot; ki++)
      //{
      //   package_imslconlin_ps->lowbounds_dv->v[ki] = 1.0E-10; //Bounded away from 0.0, because Dirichlet gives -infty for logPrior.
      //   package_imslconlin_ps->upperbounds_dv->v[ki] = 1.0;
      //}
      //--- For all model free parameters (bar transition matrix) to be positive -- AUTOMATED!  If negative, one needs to overwrite these.
      //--- This loop may be overwritten by the following package_imslconlin_ps->lowbounds_dv->v.
      //--- The upper bound is already specified in CreateTSpackagae_imslconlin() in optpackage.c.
      //for (ki=dsgemodel_ps->xphi_dv->n-1; ki>=0; ki--)
      //   package_imslconlin_ps->lowbounds_dv->v[ki] = 1.0E-10; //Bounded away from 0.0 to be save because IMSL may go over 0.0
                                //which gives a bad starting value for successive runs.
      //--- Specific bounds from datainp_common.prn -- AUTOMATED!
      // if (dsgemodel_ps->indx_hard_restrict==0) //hard0_ tag: default model (persistence for all exogenous shocks)
      // {
      //    if ((dsgemodel_ps->hard0_simple_locs_iv->n != dsgemodel_ps->hard0_simple_lowvals_dv->n) || (dsgemodel_ps->hard0_simple_locs_iv->n != dsgemodel_ps->hard0_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard0_simple_locs_iv, hard0_simple_lowvals_dv,"
      //                       "\t and hard0_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard0_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard0_simple_locs_iv->v[ki]] = dsgemodel_ps->hard0_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard0_simple_locs_iv->v[ki]] = dsgemodel_ps->hard0_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==1) //hard1_ tag: pre-fixed glambdastar and geta.
      // {
      //    if ((dsgemodel_ps->hard1_simple_locs_iv->n != dsgemodel_ps->hard1_simple_lowvals_dv->n) || (dsgemodel_ps->hard1_simple_locs_iv->n != dsgemodel_ps->hard1_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard1_simple_locs_iv, hard1_simple_lowvals_dv,"
      //                       "\t and hard1_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard1_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard1_simple_locs_iv->v[ki]] = dsgemodel_ps->hard1_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard1_simple_locs_iv->v[ki]] = dsgemodel_ps->hard1_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==2) //hard2_ tag: no persistence on all exogenous processes.
      // {
      //    if ((dsgemodel_ps->hard2_simple_locs_iv->n != dsgemodel_ps->hard2_simple_lowvals_dv->n) || (dsgemodel_ps->hard2_simple_locs_iv->n != dsgemodel_ps->hard2_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard2_simple_locs_iv, hard2_simple_lowvals_dv,"
      //                       "\t and hard2_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard2_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard2_simple_locs_iv->v[ki]] = dsgemodel_ps->hard2_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard2_simple_locs_iv->v[ki]] = dsgemodel_ps->hard2_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==3) //hard3_ tag:
      // {
      //    if ((dsgemodel_ps->hard3_simple_locs_iv->n != dsgemodel_ps->hard3_simple_lowvals_dv->n) || (dsgemodel_ps->hard3_simple_locs_iv->n != dsgemodel_ps->hard3_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard3_simple_locs_iv, hard3_simple_lowvals_dv,"
      //                       "\t and hard3_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard3_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard3_simple_locs_iv->v[ki]] = dsgemodel_ps->hard3_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard3_simple_locs_iv->v[ki]] = dsgemodel_ps->hard3_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==103) //hard103_ tag:
      // {
      //    if ((dsgemodel_ps->hard103_simple_locs_iv->n != dsgemodel_ps->hard103_simple_lowvals_dv->n) || (dsgemodel_ps->hard103_simple_locs_iv->n != dsgemodel_ps->hard103_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard103_simple_locs_iv, hard103_simple_lowvals_dv,"
      //                       "\t and hard103_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard103_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard103_simple_locs_iv->v[ki]] = dsgemodel_ps->hard103_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard103_simple_locs_iv->v[ki]] = dsgemodel_ps->hard103_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==203) //hard203_ tag:
      // {
      //    if ((dsgemodel_ps->hard203_simple_locs_iv->n != dsgemodel_ps->hard203_simple_lowvals_dv->n) || (dsgemodel_ps->hard203_simple_locs_iv->n != dsgemodel_ps->hard203_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard203_simple_locs_iv, hard203_simple_lowvals_dv,"
      //                       "\t and hard203_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard203_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard203_simple_locs_iv->v[ki]] = dsgemodel_ps->hard203_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard203_simple_locs_iv->v[ki]] = dsgemodel_ps->hard203_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==4) //hard4_ tag: no persistence on all exogenous processes.
      // {
      //    if ((dsgemodel_ps->hard4_simple_locs_iv->n != dsgemodel_ps->hard4_simple_lowvals_dv->n) || (dsgemodel_ps->hard4_simple_locs_iv->n != dsgemodel_ps->hard4_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard4_simple_locs_iv, hard4_simple_lowvals_dv,"
      //                       "\t and hard4_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard4_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard4_simple_locs_iv->v[ki]] = dsgemodel_ps->hard4_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard4_simple_locs_iv->v[ki]] = dsgemodel_ps->hard4_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==5) //hard5_ tag:
      // {
      //    if ((dsgemodel_ps->hard5_simple_locs_iv->n != dsgemodel_ps->hard5_simple_lowvals_dv->n) || (dsgemodel_ps->hard5_simple_locs_iv->n != dsgemodel_ps->hard5_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard5_simple_locs_iv, hard5_simple_lowvals_dv,"
      //                       "\t and hard5_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard5_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard5_simple_locs_iv->v[ki]] = dsgemodel_ps->hard5_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard5_simple_locs_iv->v[ki]] = dsgemodel_ps->hard5_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==105) //hard105_ tag:
      // {
      //    if ((dsgemodel_ps->hard105_simple_locs_iv->n != dsgemodel_ps->hard105_simple_lowvals_dv->n) || (dsgemodel_ps->hard105_simple_locs_iv->n != dsgemodel_ps->hard105_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard105_simple_locs_iv, hard105_simple_lowvals_dv,"
      //                       "\t and hard105_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard105_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard105_simple_locs_iv->v[ki]] = dsgemodel_ps->hard105_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard105_simple_locs_iv->v[ki]] = dsgemodel_ps->hard105_simple_highvals_dv->v[ki];
      //    }
      // }
      // else if (dsgemodel_ps->indx_hard_restrict==205) //hard205_ tag:
      // {
      //    if ((dsgemodel_ps->hard205_simple_locs_iv->n != dsgemodel_ps->hard205_simple_lowvals_dv->n) || (dsgemodel_ps->hard205_simple_locs_iv->n != dsgemodel_ps->hard205_simple_highvals_dv->n))
      //       fn_DisplayError("main.c: dimensions of hard205_simple_locs_iv, hard205_simple_lowvals_dv,"
      //                       "\t and hard205_simple_highvals_dv must match");
      //    for (ki=dsgemodel_ps->hard205_simple_locs_iv->n-1; ki>=0; ki--)
      //    {
      //       package_imslconlin_ps->lowbounds_dv->v[dsgemodel_ps->hard205_simple_locs_iv->v[ki]] = dsgemodel_ps->hard205_simple_lowvals_dv->v[ki];
      //       package_imslconlin_ps->upperbounds_dv->v[dsgemodel_ps->hard205_simple_locs_iv->v[ki]] = dsgemodel_ps->hard205_simple_highvals_dv->v[ki];
      //    }
      // }
      // else  fn_DisplayError("main():  No other cases specified for IMSL simple bounds");
   }
   #endif


   //-----------------
   // Main matter.
   //-----------------
   time(&dsgemodel_ps->prog_begtime);  //Beginning time of the whole program.
   InitializeGlobleSeed(dsgemodel_ps->randomseed);
   csminwel_randomseedChanged(dsgemodel_ps->randomseed);  //Using the same (or different) seednumber to reset csminwel seednumber for random perturbation.
   initialize_generator(dsgemodel_ps->randomseed); //For my own random number generator in rand.c. Used by, say, betarnd() in Draw_xphi_dvFromPrior_Specific().
   //--- Finding the log scale factor of the prior density that has zeros in the restricted parameter region.
   if (dsgemodel_ps->scale4logpriordensity<0.0)  //If <0.0, means that the scale needs to be computed.  Otherwise, it will be supplied in datainp_TAGE.prn.
      dsgemodel_ps->scale4logpriordensity = tz_GetScaleForLogpriordensity(smodel_ps);
   //=== Finding the peak value of the logLH or logPosterior
   if (dsgemodel_ps->indxEstFinder)
   {
      //======= Minimization problem =======
      for (drawi=dsgemodel_ps->nDrawsInitsp-1; drawi>=0; drawi--)
      {
         //--- Step 2.
         InitializeForMinproblem(n_args_cl, args_cl, minpack_ps, filename_inpdatasp_cv->v); //Initialization for minimization.


         //--------------------------
         //-- Step 3: npsolve first and then csminwel.
         //--  csminwel minimization where
         //--  minpack_ps->x_dv contains the minimized vector of parameters.
         //--  minpack_ps->fret contains the minimized value.
         //--------------------------
         #ifdef _NPSOL_   
         //=== NPSOL kicks in.
         if (dsgemodel_ps->Use_NPSOL_opt)
         {
           struct TSpackage_npsolconlin_tag *package_npsolconlin_ps = CreateTSpackagae_npsolconlin(dsgemodel_ps->nfreepars_tot, dsgemodel_ps->npsol_nlin, dsgemodel_ps->npsol_nclin);

           // TODO: deal with linear/nonlinear constraints

           minfinder_noblocknpsolconlin(package_npsolconlin_ps, minpack_ps, filename_sp_vec_minproj, dsgemodel_ps->xphi_dv->n);
           CopyVector0(minpack_ps->x_dv, package_npsolconlin_ps->xsaved_dv);
         }
         #endif

         //=== csminwel kicks in.   
         if (dsgemodel_ps->indx_const)  minfinder(minpack_ps);  //Constant-parameter case.
         else  minfinder_blockcsminwel(minpack_ps, indxFindMLE);  //Time-varying case.

         #if defined(_IMSLOPTPACKAGE_)
         if (dsgemodel_ps->Use_IMSL_opt)
         {
            //--------------------------
            //-- Step 3: IMSL linearly-constrained minimization where
            //--  minpack_ps->x_dv contains the minimized vector of parameters
            //--  minpack_ps->fret can be obtained by calling csminwell again or LogLikelihood_StatesIntegratedOut().
            //--------------------------
            minfinder_noblockimslconlin(package_imslconlin_ps, minpack_ps, filename_sp_vec_minproj, dsgemodel_ps->xphi_dv->n);
            CopyVector0(minpack_ps->x_dv, package_imslconlin_ps->xsaved_dv);

            //--------------------------
            //-- Step 3: csminwel minimization AGAIN where
            //--  minpack_ps->x_dv contains the minimized vector of parameters.
            //--  minpack_ps->fret contains the minimized value.
            //--------------------------
            if (dsgemodel_ps->indx_const)  minfinder(minpack_ps);  //Constant-parameter case.
            else  minfinder_blockcsminwel(minpack_ps, indxFindMLE);  //Time-varying case.
         }
         #endif

         //======= Writing to a file log posterior kernel and estimated results =======
         //--- Refreshing the specific-model estimates and final individual parameters from minpack_ps->x_dv.
         if (dsgemodel_ps->indx_const)  logOverallPosteriorKernal_const(smodel_ps, minpack_ps->x_dv);  //Constant-parameter.
         else
         {
            //--- The following two functions are DW's.
            //SetupObjectiveFunction(smodel_ps, minpack_ps->x_dv->v, minpack_ps->x_dv->v+dsgemodel_ps->xphi_dv->n, minpack_ps->x_dv->v);
            //logval = PosteriorObjectiveFunction(minpack_ps->x_dv->v, minpack_ps->x_dv->n);  //Refreshing.
            //--- Or better alternatives.
            ConvertFreeParametersToTheta(smodel_ps, minpack_ps->x_dv->v);
            ConvertFreeParametersToQ(smodel_ps, minpack_ps->x_dv->v+NumberFreeParametersTheta(smodel_ps));
            ThetaChanged(smodel_ps); //DW function, which triggers to reset the forward recursion at the beging and set my refreshed tracker to 0.
            dsgemodel_ps->LHvalatpostpeak = LogLikelihood_StatesIntegratedOut(smodel_ps);  //log(P(Y_T|theta, Q)).  Waggoner's function.
            dsgemodel_ps->peaklogpost = dsgemodel_ps->LHvalatpostpeak + LogPrior(smodel_ps);  //logPosterior: DW's function.
         }
         //--- Writing log posterior kernel and estimated results for each starting point to a row of the file outdataest_localmodes_TAG.prn.
         ftd_FprintEstResults(fptr_kernel7x_dv_localmodes, minpack_ps);
      }   
   }
   else  InitializeForMinproblem(n_args_cl, args_cl, minpack_ps, filename_inpdatasp_cv->v);




   //====================== The following order matters! =================================//
   //--- [1] Refreshing the specific-model estimates and final individual parameters from minpack_ps->x_dv.
   if (dsgemodel_ps->indx_const)  logOverallPosteriorKernal_const(smodel_ps, minpack_ps->x_dv);  //Constant-parameter.
   else
   {
      //--- The following two functions are DW's.
      //SetupObjectiveFunction(smodel_ps, minpack_ps->x_dv->v, minpack_ps->x_dv->v+dsgemodel_ps->xphi_dv->n, minpack_ps->x_dv->v);
      //logval = PosteriorObjectiveFunction(minpack_ps->x_dv->v, minpack_ps->x_dv->n);  //Refreshing.
      //--- Or better alternatives.
      ConvertFreeParametersToTheta(smodel_ps, minpack_ps->x_dv->v);
      ConvertFreeParametersToQ(smodel_ps, minpack_ps->x_dv->v+NumberFreeParametersTheta(smodel_ps));
      ThetaChanged(smodel_ps); //DW function, which triggers to reset the forward recursion at the beging and set my refreshed tracker to 0.
      dsgemodel_ps->LHvalatpostpeak = LogLikelihood_StatesIntegratedOut(smodel_ps);  //log(P(Y_T|theta, Q)).  Waggoner's function.
      dsgemodel_ps->peaklogpost = dsgemodel_ps->LHvalatpostpeak + LogPrior(smodel_ps);  //logPosterior: DW's function.
   }

   //--- [2] Marginal posterior probabilities of states (i.e., given Y_T).
   if (!dsgemodel_ps->indx_const)  //Regime-switching case
   {
      postprobs_sdv.flag = postprobs2_sdv.flag = V_DEF;
      postprobs_sdv.n = postprobs2_sdv.n = dsgemodel_ps->PostProbS_dm->nrows;

      //--- Getting probs of grand regimes.
      dsgemodel_ps->dw_ProbStates_dv->flag = V_DEF;
      for (si=smodel_ps->sv->nstates-1; si>=0; si--)
      {

         postprobs_sdv.v = dsgemodel_ps->PostProbS_dm->M + mos(0, si, postprobs_sdv.n);
         if (!ProbabilitiesState(dsgemodel_ps->dw_ProbStates_dv, si, smodel_ps))
            fn_DisplayError("main.c: error (probably related to dimension) occurs when calling DW's ProbabilitiesState()");
         CopySubvector(&postprobs_sdv, 0, dsgemodel_ps->dw_ProbStates_dv, dsgemodel_ps->nlags, dsgemodel_ps->fss+1);
      }
      dsgemodel_ps->PostProbS_dm->flag = M_GE;

      //--- Getting probs of base regimes.
      InitializeConstantCell_lf(dsgemodel_ps->BaseProbS_dc, 0.0);
      for (svi=smodel_ps->sv->n_state_variables-1; svi>=0; svi--)
      {
         for (basei=smodel_ps->sv->state_variable[svi]->nbasestates-1; basei>=0; basei--)
         {
            postprobs2_sdv.v = dsgemodel_ps->BaseProbS_dc->C[svi]->M + mos(0, basei, postprobs_sdv.n);
            for (si=smodel_ps->sv->nstates-1; si>=0; si--)
            {
               si_subgrand = smodel_ps->sv->index[si][svi];
               if (smodel_ps->sv->state_variable[svi]->nlags_encoded >= 1)
                  s_t_i = smodel_ps->sv->state_variable[0]->lag_index[si_subgrand][0]; //[0] means s_t and [1] means s_{t-1}.
               else
                  s_t_i = si_subgrand;

               if (s_t_i==basei)
               {
                  postprobs_sdv.v = dsgemodel_ps->PostProbS_dm->M + mos(0, si, postprobs_sdv.n);
                  VectorPlusMinusVectorUpdate(&postprobs2_sdv, &postprobs_sdv, 1.0);
               }

            }
         }
         dsgemodel_ps->BaseProbS_dc->C[svi]->flag = M_GE;
      }
   }


   //=== Other useful computations.
   //--- After [1], computing the covariance matrix or Hessian for standard errors at the log posterior peak, using the outer-product Hessian.
   if (dsgemodel_ps->flag_ComputingHessian)
   {
      //--- Inner product seems not to give us a good result.
      //ComputeHessianFrom2ndDerivative(dsgemodel_ps->Hessian_dm, smodel_ps, minpack_ps->x_dv);
      //--- Outer product.
      ComputeHessianFromOuterProduct(dsgemodel_ps->Hessian_dm, smodel_ps, minpack_ps->x_dv);
      //--- We don't compute the following becaue Hessian may NOT be SPD, and thus chol() will fail.
      //ComputeCovarianceFromOuterProduct(dsgemodel_ps->Omega_dm, smodel_ps, minpack_ps->x_dv);
   }


   //--- Ending time of the whole program.
   time(&dsgemodel_ps->prog_endtime);

   //=== Writes out the final output.
   ftd_FprintOutput(fptr_output, minpack_ps);
   //+
   ftd_FprintMatlab(fptr_matlab, smodel_ps); //Probabilities of regimes.
   ftd_FprintMatlab1(fptr_matlab1, smodel_ps); //Conditional likelihoods and forecasting errors.
   ftd_FprintMatlab2(fptr_matlab2, smodel_ps);  //For debugging purpose -- printing out the likelihood at time t and regime st.
   //ftd_const_CounterFactual(fptr_matlab3, smodel_ps); //Computing structural shocks conditional on one particular regime throughout.  ?????? Questionable historical decomps based on Junior's notes.


   printf("\n\n====================================================\n");
   printf("      Done with Tao's program!");
   printf("\n====================================================\n\n");
   fflush(stdout);
   
   //************************************************************
   //* DW simulations incluidng historical decomps and MCMC.
   //************************************************************
   // -forecast
   // -ir             (impulse responses)
   // -historical     (historical decompositions)
   // -smoothedshocks
   // -smoothedstates
   TStateModel *state_space_model = CreateStateModel_MSStateSpace(smodel_ps,1);

   /* Debugging */
   TVector parms=CreateVector(NumberFreeParametersTheta(state_space_model)+NumberFreeParametersQ(state_space_model));

   printf("**************** The following Tao and Dan implementation printouts are for debugging purposes.  *************\n");
   printf("Tao implementation - log posterior initial: %le\n",LogPosterior_StatesIntegratedOut(smodel_ps));
   printf("Tao implementation - log likelihood initial: %le\n",LogLikelihood_StatesIntegratedOut(smodel_ps));
   printf("Tao implementation - log prior theta initial: %le\n",LogPrior_theta(smodel_ps));
   printf("Tao implementation - log prior q initial: %le\n",LogPrior_q(smodel_ps));
   printf("Tao implementation - parameters:\n");
   ConvertThetaToFreeParameters(smodel_ps,pElementV(parms));
   ConvertQToFreeParameters(smodel_ps,pElementV(parms)+NumberFreeParametersTheta(smodel_ps));
   dw_PrintVector(stdout,parms,"%lf ");

   printf("Dan implementation - log posterior initial: %le\n",LogPosterior_StatesIntegratedOut(state_space_model));
   printf("Dan implementation - log likelihood initial: %le\n",LogLikelihood_StatesIntegratedOut(state_space_model));
   printf("Dan implementation - log prior theta initial: %le\n",LogPrior_theta(state_space_model));
   printf("Dan implementation - log prior q initial: %le\n",LogPrior_q(state_space_model));
   printf("Dan implementation - parameters:\n");
   ConvertThetaToFreeParameters(state_space_model,pElementV(parms));
   ConvertQToFreeParameters(state_space_model,pElementV(parms)+NumberFreeParametersTheta(state_space_model));
   dw_PrintVector(stdout,parms,"%lf ");

   /* // debug stuff */
   /* T_MSStateSpace *statespace=(T_MSStateSpace*)(state_space_model->theta); */
   /* int ii; */
   /* for (ii=0; ii < statespace->nbasestates; ii++) */
   /*   { */
   /*     printf("regime %d\n",ii); */
   /*     dw_PrintVector(stdout,statespace->b[ii],"%le "); */
   /*   } */

   /**/

   /*--- DW's MCMC program.
   dw_command_line_statespace_output(n_args_cl,args_cl,state_space_model);  //Does everything: historical decompps and MCMC.
   */

   printf("\n\nDone with Dan's program!\n\n");
   fflush(stdout);

   // Equi-Energy sampling starts from here
	cout << "equil energy sampling started: \n"; 
	//=== Parameters used for equi-energy sampling === //

	// Initialize MPI //
	MPI_Init(&n_args_cl, &args_cl);
        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Equi-Energy Model
	CEquiEnergy_TState model;
	//target model; 
	model.target_model = state_space_model;	
	model.SaveTargetModelOriginalSetting(); 
	model.parameter = new CEESParameter; 
	
	// Equi-energy sampling parameters
	// directory to store samples
        model.parameter->storage_dir = getenv("HOME")+string("/equal_energy_generic_interface/result/"); 
        model.parameter->storage_marker = 10000;
	stringstream run_id_string; 
	run_id_string.str(string()); 
	run_id_string << time(NULL); 
        model.parameter->run_id = string(dw_ParseString_String(n_args_cl, args_cl, "RunID", run_id_string.str().c_str()));


	// highest and lowest levels
	int if_original = dw_ParseInteger_String(n_args_cl, args_cl, "Original", 0);
	if (if_original == 0)	// use multiple energy levels
	{
		model.parameter->number_energy_level = dw_ParseInteger_String(n_args_cl, args_cl, "nLevel", 10);
        	model.parameter->pee = 0.3;
		model.parameter->highest_level = dw_ParseInteger_String(n_args_cl, args_cl, "lHigh", model.parameter->number_energy_level-1);
       		model.parameter->lowest_level = dw_ParseInteger_String(n_args_cl, args_cl, "lLow", 0);
	}
	else // only consider one level -- the original model
	{
		model.parameter->number_energy_level = 1; 
		model.parameter->highest_level = model.parameter->lowest_level = 0; 
	}
       	model.parameter->t0 = dw_ParseFloating_String(n_args_cl, args_cl, "T0", 1.0);
	model.parameter->tk_1 = dw_ParseFloating_String(n_args_cl, args_cl, "TK", 1000); 
       	model.parameter->SetTemperature();

	// burn-in length, simulation length, and deposit frequency
	model.parameter->thin = dw_ParseInteger_String(n_args_cl, args_cl, "thin", 1);
	model.parameter->THIN = dw_ParseInteger_String(n_args_cl, args_cl, "THIN", 100); 
        model.parameter->simulation_length = dw_ParseInteger_String(n_args_cl, args_cl, "ndraws", 10000);
	model.parameter->simulation_length = model.parameter->simulation_length; 
	
	// Block size for random blocking
	// 0: no random blocking
	model.parameter->size_per_block = dw_ParseInteger_String(n_args_cl, args_cl, "BlockSize", 0); 

	if (if_original > 0)	// if only the original model is used, 
		model.if_bounded = false; 

	// model.current_sample; 
	char *mode_file = dw_ParseString_String(n_args_cl, args_cl, "ModeFile", (char *)NULL); 
	if ( !(mode_file != NULL && model.InitializeFromFile(string(mode_file))) && !(mode_file == NULL && model.InitializeFromTarget()) )
	{
		cerr << "Initialization failed.\n";  
		abort(); 
	}
	CSampleIDWeight mode = model.current_sample; 
	
	// storage 
	model.storage = new CStorageHead(mode.GetSize_Data(), my_rank, model.parameter->run_id, model.parameter->storage_marker, model.parameter->storage_dir, model.parameter->number_energy_level); 
	
	// metropolis 
	model.metropolis = new CMetropolis(&model); 
	
	// master dispatches while slave runs tasks
	dw_initialize_generator(time(NULL)+my_rank*1000);
	if (my_rank == 0)
		master_deploying(n_args_cl, args_cl, model, mode); 
	else 
		slave_computing(n_args_cl, args_cl, model, mode);

	// end
	model.RecoverTargetModelOriginalSetting(); 
        MPI_Finalize();
        cout << "equi energy sampling : simulation done.\n";	

   // Equi-Energy sampleing ends here

   //************************************************************
   //*** DW's OLD code starts from here.
   //************************************************************
   //if (dsgemodel_ps->indxSimulateMCMC)
   //{
      //For some calculations such as historical decomposition, we need to call
      //For fast (almost half of the MCMC time) computations, we should always call
      //    dw_model=CreateStateModel_MSStateSpace(smodel_ps,0);
      //In some situations (mainly hisotrical decompositions), an additional lag in the Markov process 
      //is needed to compute historical decompositions. In this case, we first call
      //    FreeStateModel(dw_model);
      //Then call
      //    dw_model=CreateStateModel_MSStateSpace(smodel_ps,1);
      //where means an additional lag, which would lead to too much of MCMC time and should be avoided in MCMC.
      //--- Create new data structure (introducing an extra lag to compute historical decomposition for switching models)
      //dw_model=CreateStateModel_MSStateSpace(smodel_ps,0);

      //Simulation command-line parameters:
      //-ndraws : number of posterior draws to save (default = 1000)
      //-burnin : number of posterior burn-in draws (default = 0.1 * ndraws)
      //-thin : thinning factor.  Total number of draws made is thin*ndraws + burnin (default = 1)
      //-thin : Thinning factor.  Only 1/thin of the draws in posterior draws file are used. (default = 1)
      //-states : Print impluse responses of state variables to shocks (default = off)
      //-observables : Print impulse responses of observable variables to shocks (default = off)
      //-ergodic : Print impuluse responses using ergodic distribution as initial condition for shocks.
      //-regimes : Print impulse responses as if each regime were permanent.
      //-regime s : Print impulse responses as if regime s were permanent.
      //--- Generates posterior draws by passing either dw_model or smodel_ps - unclear which is more efficient
      //--- It takes a lot of time to adpatively finding the optimal scale for Metropolis.
      //dw_Simulate_command_line(n_arg,args_cl,smodel_ps,minpack_ps->x_dv,0);  //Old code but works for the time being.
      //dw_Simulate_command_line(n_arg,args_cl,dw_model,minpack_ps->x_dv,0);  //New code that seems to work.
      //printf("Done with simulationg posterior draws.\n\n");


      //-percentiles n p_1 p_2 ... p_n : Percentiles to compute. The first parameter
      //   after percentiles must be the number of percentiles and the following 
      //   values are the actual percentiles. (default = 3  0.16  0.50  0.84)
      //-horizon : The horizon over which to compute the impulse responses. (default = 12)
      //--- Generates impulse responses using the posterior draws generated by dw_Simulate_command_line().
      //dw_state_space_impulse_response_command_line(n_arg,args_cl,dw_model);
      //printf("Done with computing impulse responses.\n\n");

      //-pt : proposal type (1=gaussian, 2=power (default), 3=truncated power, 4=step, 5=truncated gaussian)
      //-d  : number proposal draws (default = 100,000)
      //--- Generates MDD by passing either dw_model or smodel_ps - unclear which is more efficient
      //--- Takes time when -d (number of proposal draws) is large.
      //dw_marginal_data_density_command_line(n_arg,args_cl,dw_model,0);
      //printf("Done with computing marginal data density.\n\n");

      //-horizon : The horizon over which to compute the impulse responses. 
      //   (default = 8)
      //-error_bands : Output error bands.  (default = off - only median is computed)
      //-percentiles n p_1 p_2 ... p_n : Percentiles to compute. The first parameter
      //   after percentiles must be the number of percentiles and the following 
      //   values are the actual percentiles. (default = 3  0.16  0.50  0.84 if 
      //   error_bands flag is set and default = 1 0.50 otherwise)
      //-parameter_uncertainty : Apply parameter uncertainty when computing error 
      //   bands.
      //-shocks_per_parameter : Number of shocks and regime paths to draw for each
      //   parameter draw.  (default = 1 if parameter_uncertainty is set and 
      //   default = 10,000 otherwise)
      //-thin : Thinning factor.  Only 1/thin of the draws in posterior draws file 
      //   are used. (default = 1)
      //-regimes : Produces forecasts as if each regime were permanent.
      //-regime s : Produces forecasts as if regime s were permanent.
      //-mean : Produces mean forecast.  (default = off)
      //--- Generats forecasts with options for error bands (new forecasting code)
      //dw_state_space_forecast_command_line(n_arg,args_cl,dw_model);
      //printf("Done with forecasting.\n\n");


      //=== Old code for producing historical decompositions, smoothed shocks, and forecasts (unconditional or conditional on shocks)
      //    The following function call is in dw_temp_output.c in the subdirectory switching.
      //FreeStateModel(dw_model);
      //dw_model=CreateStateModel_MSStateSpace(smodel_ps,1);
      
      //printf("Transition output stuff that may or may not useful\n   -- for now historical decompositions may be useful\n");
      //dw_DansOutputSomeStuff(smodel_ps, cl_modeltag, minpack_ps->x_dv);

   //   printf("Done with the DW's MCMC program.\n");
   //}



   //=== Memory destruction.
   //--- The order matters.
   DestroyTSargs_blockcsminwel(args_blockcsminwel_ps);
   DestroyTSetc_minproj(etc_minproj_ps);
   DestroyTSminpack(minpack_ps);
   //+
   #if defined(_IMSLOPTPACKAGE_)
   DestroyTSpackagae_imslconlin(package_imslconlin_ps);
   #endif
   //+ To be destroyed last.
   DestroyTSdsgemodel(dsgemodel_ps);
   FreeStateModel(smodel_ps);   //Frees everything related to Waggoner's Markov-switching model.
   //
   DestroyVector_c(filename_inpdatasp_cv);
   //+
   tzFclose(fptr_rawdata);
   tzFclose(fptr_common);
   tzFclose(fptr_specific);
   tzFclose(fptr_markov);
   //+
   tzCloseDebugFile();
   tzCloseOptIntermFile();
   tzFclose(fptr_matlab);
   tzFclose(fptr_matlab1);
   tzFclose(fptr_matlab2);
   tzFclose(fptr_matlab3);
   tzFclose(fptr_kernel7x_dv_localmodes);
   tzFclose(fptr_output);


   //=== Done!
   printf("\n*****************************************");
   printf("\n*** Program is successfully executed! ***\n");
   printf("*****************************************\n");
   return (EXIT_SUCCESS);
}



//-----------------------
// Minimization problem.
//-----------------------
//------- Step 1: Only for the constant-parameter case. -------
static double minobj(struct TSminpack_tag *minpack_ps)
{
   //Returns the negative value of log PosteriorFunction.
   return ( -logOverallPosteriorKernal_const(((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->smodel_ps, minpack_ps->xtemp_dv) );
}


//------- Step 2. -------
static void InitializeForMinproblem(int n_args_cl, char **args_cl, struct TSminpack_tag *minpack_ps, char *filename_inpdatasp)
{
   //$$$$$$$ This function is made automatic (I hope) to take care of 
   //$$$$$$$    (1) constant-parameter case without using DW's functions;
   //$$$$$$$    (2) allowing us to generate parameters randomly, which depends on the specific model.
   //$$$$$$$    If the user wants to change certain features, the only changes one may make are
   //                (1) the name of the user's model structure.  Here it is TSdsgemodel_tag.
   //                (2) the name of xphi_dv;
   //                (3) the name of indxStartValuesForMin;
   //                (4) the name of indx_const;
   //                (5) the functoin name of logOverallPosteriorKernal_const();
   //$$$$$$$    Otherwise, everything else is automatic.
   //
   //Outputs:
   //  minpack_ps->x_dv and minpack_ps->xtemp_dv:
   //    The 1st xphi_dv->n elements of x_dv are model parameters (excluding those in the transition matrices).
   //    The 2nd-part or rest of the elements of x_dv are the free parameters in the transition matrices.
   //---
   //For indxStartValuesForMin (corresponding to /c in runprog.bat):
   //  0: read in from the last estimated results contained in inpdatasp_min_TAG.prn (= a row vector of free parameters only).
   //  1: starts from the fixed values for xphi_dv, read in from datainp_*.prn.
   //  2: randomly selects the initial starting values drawn from the prior density.
   //  3: randomly selects the initial starting values drawn from the simulation_TAG.out generated by DW MCMC code (each row = [logPosterior, logLH, free parameters]).

   FILE *fptr_sp_input = NULL;     //File pointer to input file containing the starting point.
   int _n, _i, nqs;
   int flag_success, nTrys, nTermination=500;
   TSdvector xphi_sdv, xqs_sdv;
   TSdvector *x_dv = minpack_ps->x_dv;
   TSdvector *x0_dv = minpack_ps->x0_dv;
   //---
   struct TStateModel_tag *smodel_ps = (struct TStateModel_tag *)((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->smodel_ps;
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   int indxStartValuesForMin = dsgemodel_ps->indxStartValuesForMin;
   int indx_const = dsgemodel_ps->indx_const;
   int nphis = dsgemodel_ps->xphi_dv->n;
   //--- For Markov-switching model.
   int n1, n2;
   double *x1_pd, *x2_pd;
   //--- For option /c 3 only (using DW's code)
   char *filename_mcmc=(char *)NULL;  //Memory will be not created here because it'll simply point to the filename from the command line.  Thus, no need to destroy this.  
   int row, n_theta=NumberFreeParametersTheta(smodel_ps), n_q=NumberFreeParametersQ(smodel_ps);
   PRECISION logposterior, loglikelihood;
   TMatrix draws=(TMatrix)NULL;
   TVector parameters=(TVector)NULL;
   
   
   //=== Odds and ends.
   TSdvector *tfxqs_dv = NULL;

   xphi_sdv.flag = V_DEF;
   xphi_sdv.n = nphis;
   xphi_sdv.v = x_dv->v;
   if (!indx_const)   //Regime switching case.
   {
      if ( (nqs=NumberFreeParametersQ(smodel_ps)) != (x_dv->n - nphis) )
         fn_DisplayError(".../*_estmcmc.c/InitializeForMinproblem(): Minimization vector must have length\n"
                         "     equal to # of free model parameters plus # of free transition matrix parameters");

      xqs_sdv.flag = V_DEF;
      xqs_sdv.n = nqs;
      xqs_sdv.v = x_dv->v + xphi_sdv.n;
      //
      tfxqs_dv = CreateVector_lf(xqs_sdv.n);
   }

   if (indxStartValuesForMin == 1)
   {
      //=== Copying xphi_dv to the minimization routine.
      if (!indx_const)  //Regime-switching case.
      {
         CopyVector0(&xphi_sdv, dsgemodel_ps->xphi_dv);
         if (1)   //Use Q (transitiona matrices) in the minimization problem.
            ConvertQToFreeParameters(smodel_ps, xqs_sdv.v);     //Waggnoer's function.
         else   //Use log(Q) (log(transitiona matrices)) in the minimization problem.
         {
            ConvertQToFreeParameters(smodel_ps, tfxqs_dv->v);     //Waggnoer's function.
            tfxqs_dv->flag = V_DEF;
            ConvertVector2log(&xqs_sdv, tfxqs_dv);   //Converting it from normal value to log value for minimization.
         }
         x_dv->flag = V_DEF;
      }
      else  CopyVector0(x_dv, dsgemodel_ps->xphi_dv);   //Constant-parameter case.
   }
   else if (!indxStartValuesForMin)   //Reading from the last estimated values stored in the inpdatasp_min_*.prn file.
   {
      if (!(fptr_sp_input = fopen(filename_inpdatasp,"r")))
      {
         printf("\n\n*_estmcmc.c/InitializeForMinproblem(): With /c 0 in the command line, the file %s must exist for reading free parameter values!\n", filename_inpdatasp);
         fprintf(tzGetDebugFile(), "\n\n*_estmcmc.c/InitializeForMinproblem(): The file %s must exist for reading free parameter values!\n", filename_inpdatasp);
         exit(EXIT_FAILURE);
      }
      rewind(fptr_sp_input);   //To be safe by putting the pointer at the beginning of the file.

      for (_n=x_dv->n, _i=0; _i<_n; _i++)
         if (fscanf(fptr_sp_input, " %lf ", x_dv->v+_i) != 1)
         {
            printf("Error: *_estmcmc.c/InitializeForMinproblem() -- cannot read the number from the file inpdatasp_min_TAG.prn.\n");
            exit(EXIT_FAILURE);
         }
      x_dv->flag = V_DEF;

      if (!indx_const)  //Regime-switching case.
      {
         CopyVector0(dsgemodel_ps->xphi_dv, &xphi_sdv);
         ConvertFreeParametersToTheta(smodel_ps, x_dv->v);
         if (1)   //Use Q (transitiona matrices) in the minimization problem.
            ConvertFreeParametersToQ(smodel_ps, xqs_sdv.v); //This DW function always sets the flag that q has changed.
         else   //Use log(Q) (log(transitiona matrices)) in the minimization problem.
         {
            ConvertVector2exp(tfxqs_dv, &xqs_sdv);   //Converting it from log value back to the normal value.
            ConvertFreeParametersToQ(smodel_ps, tfxqs_dv->v);
         }
      }
      else
         CopyVector0(dsgemodel_ps->xphi_dv, x_dv);
   }
   else if (indxStartValuesForMin == 2)
   {
      //--- Getting a draw of xphi_dv from the prior distribution.
      flag_success = 0;
      nTrys = 0;
      while (!flag_success)
      {
         Draw_xphi_dvFromPrior_Specific(smodel_ps);
         flag_success = RefreshEverything(smodel_ps);
         if (nTrys++ > nTermination)  fn_DisplayError("*_estmcmc.c/InitializeForMinproblem() in main.c: cannot find a legitimate draw from the prior");
      }
      dsgemodel_ps->xphi_dv->flag = V_DEF;


      //=== Copying xphi_dv to the minimization routine.
      if (!indx_const)  //Regime-switching case.
      {
         CopyVector0(&xphi_sdv, dsgemodel_ps->xphi_dv);
         if (1)   //Use Q (transitiona matrices) in the minimization problem.
            ConvertQToFreeParameters(smodel_ps, xqs_sdv.v);     //Waggnoer's function.
         else   //Use log(Q) (log(transitiona matrices)) in the minimization problem.
         {
            ConvertQToFreeParameters(smodel_ps, tfxqs_dv->v);     //Waggnoer's function.
            tfxqs_dv->flag = V_DEF;
            ConvertVector2log(&xqs_sdv, tfxqs_dv);   //Converting it from normal value to log value for minimization.
         }
         x_dv->flag = V_DEF;
      }
      else  CopyVector0(x_dv, dsgemodel_ps->xphi_dv);   //Constant-parameter case.
   }
   else if (indxStartValuesForMin == 3)
   {
      //======= Random starting values selected from simulation_TAG*.out generated by DW.  Piggy-bag on DW's code.  
      //--- Initialize generator.  No need to because we use TZ's unirnd() instead of dw_uniform_rnd().
      ////dw_initialize_generator(dw_ParseInteger_String(n_args_cl,args_cl,"seed",dsgemodel_ps->randomseed));
      //--- Create parameters vector
      parameters=CreateVector(2+n_theta+n_q);  //2+ because the first 2 elements are logPosterior and logLH.
     
      if (filename_mcmc=dw_ParseString_String(n_args_cl, args_cl, "random",(char*)NULL))
         if (draws=dw_ReadPosteriorDraws((FILE*)NULL,filename_mcmc,(char*)NULL,n_theta+n_q))
         {
            //row=(int)floor(dw_uniform_rnd()*RowM(draws));
            row=(int)floor(unirnd()*RowM(draws));
            if (row >= RowM(draws)) row=RowM(draws)-1;
            RowVector(parameters,draws,row);  //Copying to ``parameters'' the drawn parameter values with logPosterior and logLH being the first 2 elements.
            
            if (n_theta) 
            {
               memcpy(xphi_sdv.v, pElementV(parameters)+2, n_theta*sizeof(double));  //Warning: PRECISION may be incompatible with double; in this case, set PRECISION to double in DW's file prcsn.h.
               ConvertFreeParametersToTheta(smodel_ps,pElementV(parameters)+2);
            }   
            CopyVector0(dsgemodel_ps->xphi_dv, &xphi_sdv);
            if (n_q) 
            {
               memcpy(xqs_sdv.v, pElementV(parameters)+2+n_theta, n_q*sizeof(double));  //Warning: PRECISION may be incompatible with double; in this case, set PRECISION to double in DW's file prcsn.h.
               ConvertFreeParametersToQ(smodel_ps, pElementV(parameters)+2+n_theta);
            }   
            x_dv->flag = V_DEF;
            printf("\nUsing row %d as random starting value from %s\n",row,filename_mcmc);
            printf("Initial log posterior = %le (recomputed) -  %le (read in)\n",logposterior=LogPosterior_StatesIntegratedOut(smodel_ps),ElementV(parameters,0));
            printf("Initial log likelihood = %le (recomputed) -  %le (read in)\n",loglikelihood=LogLikelihood_StatesIntegratedOut(smodel_ps),ElementV(parameters,1));
            //ElementV(parameters,0)=logposterior;   //Reset.
            //ElementV(parameters,1)=loglikelihood;  //Reset.
            
            //--- Destroying memory.
            FreeMatrix(draws);
         }
      else
      {
         printf("Error: *_estmcmc.c/InitializeForMinproblem() -- unable to open the MCMC file %s\n",filename_mcmc);
         exit(EXIT_FAILURE);
      }
      //--- Destroying memory.
      FreeVector(parameters);
   }
   else  fprintf(tzGetDebugFile(), "*_estmcmc.c/InitializeForMinproblem(): the case indxStartValuesForMin (%d > 3) has not been programmed yet",indxStartValuesForMin);

//   while (logpriordensity(smodel_ps) <= (-NEARINFINITY+1.0E-12))
//   {
//      Draw_xphi_dvFromPrior_Specific(smodel_ps);
//      dsgemodel_ps->nViolates++;
//      dsgemodel_ps->nDrawsFromPrior++;
//   }

   //--- Initial or starting values of the parameters.
   CopyVector0(x0_dv, x_dv);
   ThetaChanged(smodel_ps); //Conservatively setting a flag that my model theta has changed.

   if (indx_const)  minpack_ps->fret0 = minpack_ps->fret = -logOverallPosteriorKernal_const(smodel_ps, x0_dv);  //Constant-parameter case.
   else  //Regime-switching case.
   {
      //--- Setup.
      n1 = NumberFreeParametersTheta(smodel_ps);       //Number of free model parameters.
      n2 = NumberFreeParametersQ(smodel_ps);   //Number of free transition matrix elements.
      if (x_dv->n != (n1 + n2))  fn_DisplayError("*_estmcmc.c/InitializeForMinproblem(): total number of free parameters"
                 "  must be equal to number of free model parameters + number of free q's");
      x1_pd = x_dv->v;
      x2_pd = x_dv->v + n1;
      //+

      //=== For debugging only.
      // smodel_ps->ValidForwardRecursion = 0;  //Reset so recursion can do through.
      minpack_ps->fret0 = LogPosterior_StatesIntegratedOut(smodel_ps);  //log(P(Y_T|theta, Q)).  Waggoner's function.
      // minpack_ps->fret0 = LogPrior(smodel_ps);  //logPosterior: DW's function.

      SetupObjectiveFunction(smodel_ps, x1_pd, x2_pd, x_dv->v);
      minpack_ps->fret0 = minpack_ps->fret = PosteriorObjectiveFunction(x_dv->v, x_dv->n);  //Refreshing. logPosterirPdf.  DW function.
   }

//   if (minpack_ps->fret0 >= 0.5*NEARINFINITY)  //0.5* is to make is safer to terminate the program when the bad likelihood may be slightly smaller than NEARINFINITY.
//      fn_DisplayError("*_estmcmc.c/InitializeForMinproblem(): Bad initialization. All parameters must be in the reasonable range");

   DestroyVector_lf(tfxqs_dv);
   tzFclose(fptr_sp_input);
}




//-------------------
// Prints final results to output files for Matlab and for easy eye examination.
//-------------------
#define RFORMAT  " %10.5f\n "   //R: return carriage.
#define NFORMAT  " %10.5f "     //N: no return carriage.
#define AFORMAT  " %.16e "      //A: accurate (most accurate) with no carriage return.
#define ARFORMAT  " %.16e\n "      //A: accurate (most accurate) with carriage return.
//static void ftd_const_CounterFactual(FILE *fptr_matlab3, struct TStateModel_tag *smodel_ps)
//{
//   //Works only for the constant-regime case because updating Pt_tm1_d4->F[ti]->C[0] and zt_tm1 can be tricky when regimes switch.
//   //Computing structural shocks conditional on one particular regime throughout.
//   //Have NOT done counterfactual paths under different assumptions (MANUALLY done).

//   int st = 0; //$$##  Regime is fixed for the time being.
//   int nstarts = 0; //Starting point after nlags.
//   int ti, tibeg, tibase0, t0i;
//   double loglh;
//   //--- Taken out of Waggoner's model.
//   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
//   struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps = dsgemodel_ps->kalfilmsinputs_1stapp_ps;
//   //--- Accessible variables.
//   int ny = kalfilmsinputs_1stapp_ps->ny;
//   int nz = kalfilmsinputs_1stapp_ps->nz;
//   TSdvector *xphi_dv =  dsgemodel_ps->xphi_dv;
//   TSdfourth *Pt_tm1_d4 = kalfilmsinputs_1stapp_ps->Pt_tm1_d4;
//   TSdcell *Ht_dc = kalfilmsinputs_1stapp_ps->Ht_dc; //ny-by-nz-by-nst.
//   TSdmatrix *at_dm = kalfilmsinputs_1stapp_ps->at_dm;   //ny-by-nst.
//   TSdcell *Ft_dc = kalfilmsinputs_1stapp_ps->Ft_dc; //nz-by-nz-by-nst.
//   TSdcell *Vt_dc = kalfilmsinputs_1stapp_ps->Vt_dc;     //nz-by-nz-by-nst.  Covariance matrix for the state equation.
//   TSdmatrix *bt_dm = kalfilmsinputs_1stapp_ps->bt_dm;   //nz-by-nst.
//   TSdcell *Rt_dc = kalfilmsinputs_1stapp_ps->Rt_dc;     //ny-by-ny-by-nst.  Covariance matrix for the measurement equation.
//   TSdcell *Gt_dc = kalfilmsinputs_1stapp_ps->Gt_dc;     //nz-by-ny-by-nst.  Cross-covariance.
//   //
//   TSdcell *zt_tm1_dc = kalfilmsinputs_1stapp_ps->zt_tm1_dc; //nz-by-nst-by-(nSample+1).
//   TSdmatrix *yt_dm = kalfilmsinputs_1stapp_ps->yt_dm;         //ny-by-nSample.
//   TSdcell *etdata_dc = kalfilmsinputs_1stapp_ps->etdata_dc;         //ny-by-nst-by-(fss-nstarts).
//   TSdvector etdata_sdv, yt_sdv, at_sdv, bt_sdv, zt_tm1_sdv; //yt_tm1_sdv;
//   TSdvector dates_effective_sdv;  //fss-by-1 dates
//   TSdmatrix yt_sdm;
//   //--- Accesible to memory variable but to be computed.
//   TSdmatrix *strshks_dm = dsgemodel_ps->strshks_dm; //ny-by-fss.
//   TSdvector strshks_sdv;
//   TSdmatrix strshks_sdm;
//   //===
//   TSdmatrix *HtV_dm = NULL;
//   TSdmatrix *etdata_st_dm = NULL; //ny-by-(fss-nstarts).
//   //============== For forecasts. ==============
//   TSdvector *zt_tm1_dv = NULL;
//   TSdvector *zt_dv = NULL;
//   TSdvector *ztp1_dv = NULL;
//   TSdmatrix *yconfcst_dm = NULL; //ny-by-(fss-nstarts), conditional forecasts.
//   TSdvector yconfcst_sdv;
//   //+
//   TSdvector *wny_dv = CreateVector_lf(ny);
//   TSdmatrix *Wnzbyny_dm = CreateMatrix_lf(nz,ny);
//   TSdmatrix *Wnzbynz_dm = CreateMatrix_lf(nz,nz);
//   TSdmatrix *W2nzbynz_dm = CreateMatrix_lf(nz,nz);
//   TSdmatrix *PHtran_tdata_dm = CreateMatrix_lf(nz,ny);
//   TSdvector *etdata_dv = CreateVector_lf(ny);
//   TSdmatrix *Dtdata_dm = CreateMatrix_lf(ny,ny);
//   TSdmatrix *Kt_tdata0_dm = CreateMatrix_lf(nz,ny);
//   TSdmatrix *Kt_tdata_dm = CreateMatrix_lf(nz,ny);
//   //+
//   TSdvector *w_strshks_dv = CreateVector_lf(dsgemodel_ps->fss) ; //fss
//   TSdmatrix *w_strshks_dm = CreateMatrix_lf(ny,dsgemodel_ps->fss) ; //ny-by-(fss).


//   HtV_dm = CreateMatrix_lf(ny,ny);
//   etdata_st_dm = CreateMatrix(ny,dsgemodel_ps->fss-nstarts);
//   //+
//   zt_tm1_dv = CreateVector_lf(nz);
//   zt_dv = CreateVector_lf(nz);
//   ztp1_dv = CreateVector_lf(nz);
//   yconfcst_dm = CreateMatrix_lf(ny,dsgemodel_ps->fss-nstarts); //ny-by-fss, conditional forecasts.


//   if (dsgemodel_ps->indx_const)
//   {
//      if (st>0)  fn_DisplayError("main.c/ftd_CounterFactual(): st must be no greater than 0 for the constant-parameter case");
//   }

//   dates_effective_sdv.flag = V_DEF;
//   dates_effective_sdv.n = dsgemodel_ps->fss-nstarts;
//   dates_effective_sdv.v = dsgemodel_ps->dates_dv->v + (dsgemodel_ps->nlags+nstarts);
//   //
//   bt_sdv.n = nz;
//   bt_sdv.flag = V_DEF;
//   zt_tm1_sdv.n = nz;
//   zt_tm1_sdv.flag = V_DEF;
//   //
//   etdata_sdv.n = yt_sdv.n = ny;
//   etdata_sdv.flag = yt_sdv.flag = V_DEF;
//   strshks_sdv.n = ny;
//   strshks_sdv.flag = V_DEF;
//   //
//   yt_sdm.flag = M_GE;
//   yt_sdm.nrows = ny;
//   yt_sdm.ncols =  dsgemodel_ps->fss-nstarts;
//   strshks_sdm.flag = M_GE;
//   strshks_sdm.nrows = ny;
//   strshks_sdm.ncols = dsgemodel_ps->fss-nstarts;
//   //
//   at_sdv.n = ny;
//   at_sdv.flag = V_DEF;
//   yconfcst_sdv.n = ny;
//   yconfcst_sdv.flag = V_DEF;



//   //=== Incorrect -- no longer used.  Getting structural shocks.
//   // for (ti=(dsgemodel_ps->nlags+1)+nstarts; ti<=dsgemodel_ps->nSample; ti++) //Note <= nSample because DW's function takes time as base-1.
//   // {
//   //    tibase0 = ti - 1;
//   //    t0i = ti - dsgemodel_ps->nlags - 1 - nstarts;
//   //
//   //    strshks_sdv.v = strshks_dm->M + strshks_sdv.n*(t0i+nstarts);
//   //    etdata_sdv.v = etdata_dc->C[tibase0]->M + etdata_sdv.n*st;
//   //    CopySubvector2matrix(etdata_st_dm, 0, t0i, &etdata_sdv, 0, etdata_sdv.n);
//   //
//   //    CopyVector0(&strshks_sdv, &etdata_sdv); //y_t - y_{t|t-1}.  First: backing out one-step forecast errors.
//   //    if (bdivA_rgens(&strshks_sdv, &strshks_sdv, '\\', HtV_dm))  fn_DisplayError("main.c/ftd_CounterFactual(): failed to get the structural shocks");
//   //                         //Second: backing out structural shocks.
//   // }
//   // strshks_dm->flag = M_GE;


//   //============= Conditional forecasts only (ad hoc approach). =================
//   if (1) //Counterfactual by refereshing the new ``estimates.''
//   {
//      //====== 0: MP; 1: PM; 2: WM; 3: gov spending; 4: NT; 5: preference; 6: investment tech; 7: Dep;
//      if (0)
//      {
//         //CopyMatrix0(w_strshks_dm, strshks_dm); //w_strshks_dm is never used.
//         InitializeConstantMatrix_lf(strshks_dm, 0.0);
//      }
//      else
//      {
//         InitializeConstantVector_lf(w_strshks_dv, 0.0);
//         //CopySubvector2rowmatrix(strshks_dm, 0, 0, w_strshks_dv, 0, w_strshks_dv->n); //MP
//         CopySubvector2rowmatrix(strshks_dm, 1, 0, w_strshks_dv, 0, w_strshks_dv->n); //Price Markup
//         CopySubvector2rowmatrix(strshks_dm, 2, 0, w_strshks_dv, 0, w_strshks_dv->n); //Wage Markup
//         //CopySubvector2rowmatrix(strshks_dm, 3, 0, w_strshks_dv, 0, w_strshks_dv->n); //Gov spending
//         CopySubvector2rowmatrix(strshks_dm, 4, 0, w_strshks_dv, 0, w_strshks_dv->n);   //Nentral tech
//         //CopySubvector2rowmatrix(strshks_dm, 5, 0, w_strshks_dv, 0, w_strshks_dv->n);   //Preference
//         //CopySubvector2rowmatrix(strshks_dm, 6, 0, w_strshks_dv, 0, w_strshks_dv->n); //Investment tech
//         CopySubvector2rowmatrix(strshks_dm, 7, 0, w_strshks_dv, 0, w_strshks_dv->n);   //Dep
//      }

//      ThetaChanged(smodel_ps); //DW function, which triggers to reset the forward recursion at the beging and set my refreshed tracker to 0.
//      loglh = LogLikelihood_StatesIntegratedOut(smodel_ps);  //log(P(Y_T|theta, Q)).  Waggoner's function.
//   }


//   bt_sdv.v = bt_dm->M + bt_sdv.n*st;
//   //
//   strshks_sdm.M = strshks_dm->M + strshks_sdm.nrows*nstarts;  //Used for print-out in fptr_matlab3.
//   tibeg = (dsgemodel_ps->nlags + 1) + nstarts;
//   zt_tm1_sdv.v = zt_tm1_dc->C[tibeg-1]->M + zt_tm1_sdv.n*st;
//   CopyVector0(zt_tm1_dv, &zt_tm1_sdv);
//   strshks_sdv.v = strshks_dm->M + strshks_sdv.n*(nstarts);
//   MatrixTimesVector(zt_dv,  Impact_t_exact_dc->C[st], &strshks_sdv, 1.0, 0.0, 'N'); //Impact*u_t
//   VectorPlusVector(zt_dv, &zt_tm1_sdv, zt_dv); //z_t = z_{t|t-1} + Impact*u_t.
//   for (ti=tibeg; ti<=dsgemodel_ps->nSample; ti++) //Note <= nSample because DW's function takes time as base-1.
//   {
//      tibase0 = ti - 1;
//      t0i = ti - tibeg;

////       //-------- Works.
////      //--- Observables.
////      yconfcst_sdv.v = yconfcst_dm->M + yconfcst_sdv.n*t0i;
////      MatrixTimesVector(&yconfcst_sdv, Ht_dc->C[st], &zt_tm1_sdv, 1.0, 0.0, 'N'); //H*z_{t|t-1}.
////      strshks_sdv.v = strshks_dm->M + strshks_sdv.n*(t0i+nstarts);
////      MatrixTimesVector(&yconfcst_sdv, HtV_dm, &strshks_sdv, 1.0, 1.0, 'N'); //H*z_{t|t-1} + HtV*u_t
////      at_sdv.v = at_dm->M + at_sdv.n*st;
////      VectorPlusVector(&yconfcst_sdv, &at_sdv, &yconfcst_sdv); //a + H*z_{t|t-1} + HtV*u_t.

////      //--- Getting next-period variables.
////      zt_tm1_sdv.v = zt_tm1_dc->C[ti]->M + zt_tm1_sdv.n*st;


//      //--- Computing forecast observables.
//      yconfcst_sdv.v = yconfcst_dm->M + yconfcst_sdv.n*t0i;
//      MatrixTimesVector(&yconfcst_sdv, Ht_dc->C[st], zt_dv, 1.0, 0.0, 'N'); //H*z_{t}.
//      at_sdv.v = at_dm->M + at_sdv.n*st;
//      VectorPlusVector(&yconfcst_sdv, &at_sdv, &yconfcst_sdv); //a + H*z_{t}.

//      //============ Kalman filtering condition on only one regime. ===============
//      //--- Setup.
//      MatrixTimesMatrix(PHtran_tdata_dm, Pt_tm1_d4->F[tibase0]->C[st], Ht_dc->C[st], 1.0, 0.0, 'N', 'T');

//      //--- Data.
//      //- etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata;
//      at_sdv.v = at_dm->M + st*at_dm->nrows;  //grand regime at time tbase0.
//      VectorMinusVector(etdata_dv, &yconfcst_sdv, &at_sdv);
//      MatrixTimesVector(etdata_dv, Ht_dc->C[st], zt_tm1_dv, -1.0, 1.0, 'N');
//      //+   Dtdata = Htdata*PHtran_tdata + R(:,:,tdata);
//      CopyMatrix0(Dtdata_dm, Rt_dc->C[st]);
//      MatrixTimesMatrix(Dtdata_dm, Ht_dc->C[st], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
//      ScalarTimesMatrixSquare(Dtdata_dm, 0.5, Dtdata_dm, 'T', 0.5);  //Making it symmetric against some rounding errors.
//                         //This making-symmetric is very IMPORTANT; otherwise, we will get the matrix being singular message
//                         //    and eigenvalues being negative for the SPD matrix, etc.  Then the likelihood becomes either
//                         //    a bad number or a complex number.
//      Dtdata_dm->flag = Dtdata_dm->flag | M_SU | M_SL;

//      //--- Updating zt_tm1_dm and Pt_tm1_dc for the next period.
//      if (ti < dsgemodel_ps->nSample)
//      {
//         //Updating only up to t=T-1, because the values at tp1=T will NOT be used in the likelihood function.

//         //- Kt_tdata = (Ft*PHtran_tdata+G(:,:,tdata))/Dtdata;
//         CopyMatrix0(Kt_tdata0_dm, Gt_dc->C[st]);
//         MatrixTimesMatrix(Kt_tdata0_dm, Ft_dc->C[st], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
//         BdivA_rrect(Kt_tdata_dm, Kt_tdata0_dm, '/', Dtdata_dm);
//         //+ zt_tm1(:,t) = b(:,t) + Ft*zt_tm1(:,tdata) + Kt_tdata*etdata;
//         MatrixTimesVector(zt_dv, Ft_dc->C[st], zt_tm1_dv, 1.0, 0.0, 'N');
//         MatrixTimesVector(zt_dv, Kt_tdata_dm, etdata_dv, 1.0, 1.0, 'N');
//         VectorPlusMinusVectorUpdate(zt_dv, &bt_sdv, 1.0);

//         //+ Pt_tm1(:,:,t) = Ft*Ptdata*Fttran - Kt_tdata*Dtdata*Kt_tdatatran + V(:,:,t);
//         CopyMatrix0(Pt_tm1_d4->F[ti]->C[st], Vt_dc->C[st]);
//         MatrixTimesMatrix(Wnzbyny_dm, Kt_tdata_dm, Dtdata_dm, 1.0, 0.0, 'N', 'N');
//         MatrixTimesMatrix(Wnzbynz_dm, Wnzbyny_dm, Kt_tdata_dm, 1.0, 0.0, 'N', 'T');
//         MatrixPlusMinusMatrixUpdate(Pt_tm1_d4->F[ti]->C[st], Wnzbynz_dm, -1.0);
//                               //Done with all W*_dm.
//         MatrixTimesMatrix(Wnzbynz_dm, Ft_dc->C[st], Pt_tm1_d4->F[tibase0]->C[st], 1.0, 0.0, 'N', 'N');
//         MatrixTimesMatrix(W2nzbynz_dm, Wnzbynz_dm, Ft_dc->C[st], 1.0, 0.0, 'N', 'T');
//         MatrixPlusMatrixUpdate(Pt_tm1_d4->F[ti]->C[st], W2nzbynz_dm);
//                               //Done with all W*_dm.

//         //- Updating zt_tm1_dv for next period.
//         CopyVector0(zt_tm1_dv, zt_dv);

//         //- Updating zt_dv for next period.
//         strshks_sdv.v = strshks_dm->M + strshks_sdv.n*(t0i+nstarts+1);
//         MatrixTimesVector(zt_dv, Impact_t_exact_dc->C[st], &strshks_sdv, 1.0, 1.0, 'N'); //z_{t+1|t} + Impact*u_{t+1}
//      }
//   }


//   //--- 1. Dates: 1-by-(fss-nstarts).
//   WriteVector(fptr_matlab3, &dates_effective_sdv, " %.16e ");
//   //--- 2. Actual data (for checking and debugging): ny-by-(fss-nstarts).
//   fprintf(fptr_matlab3,"\n");
//   yt_sdm.M = yt_dm->M + yt_sdm.nrows*(tibeg-1);
//   WriteMatrix(fptr_matlab3, &yt_sdm, " %.16e ");
//   //--- 3. Reduced-form residuals: ny-by-(fss-nstarts).
//   fprintf(fptr_matlab3,"\n");
//   WriteMatrix(fptr_matlab3, etdata_st_dm, " %.16e ");
//   //--- 4. Structural shocks: ny-by-(fss-nstarts).
//   fprintf(fptr_matlab3,"\n");
//   WriteMatrix(fptr_matlab3, &strshks_sdm, " %.16e ");
//   //--- 5. Conditional forecasts: : ny-by-(fss-nstarts).
//   fprintf(fptr_matlab3,"\n");
//   WriteMatrix(fptr_matlab3, yconfcst_dm, " %.16e ");



//   //===
//   DestroyMatrix_lf(HtV_dm);
//   DestroyMatrix_lf(etdata_st_dm);
//   //+
//   DestroyVector_lf(zt_tm1_dv);
//   DestroyVector_lf(zt_dv);
//   DestroyVector_lf(ztp1_dv);
//   DestroyMatrix_lf(yconfcst_dm);
//   //
//   DestroyVector_lf(wny_dv);
//   DestroyMatrix_lf(Wnzbyny_dm);
//   DestroyMatrix_lf(Wnzbynz_dm);
//   DestroyMatrix_lf(W2nzbynz_dm);
//   DestroyMatrix_lf(PHtran_tdata_dm);
//   DestroyVector_lf(etdata_dv);
//   DestroyMatrix_lf(Dtdata_dm);
//   DestroyMatrix_lf(Kt_tdata0_dm);
//   DestroyMatrix_lf(Kt_tdata_dm);
//   //
//   DestroyVector_lf(w_strshks_dv);
//   DestroyMatrix_lf(w_strshks_dm);

//}
//---
static void ftd_FprintMatlab2(FILE *fptr_matlab2, struct TStateModel_tag *smodel_ps)
{
   //For debugging purpose.
   //Printing out the likelihood at time t and regime st.

   int ti, tbase0i, si;
   //--- Taken out of Waggoner's model.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps = dsgemodel_ps->kalfilmsinputs_1stapp_ps;
   //--- Accessible variables.
   TSdvector dates_effective_sdv;  //fss-by-1 dates
   //===
   TSdvector *lht_dv = NULL; //1-by-fss.
   TSdmatrix *lhtst_dm = NULL; //h-by-fss.

   dates_effective_sdv.flag = V_DEF;
   dates_effective_sdv.n = dsgemodel_ps->fss;
   dates_effective_sdv.v = dsgemodel_ps->dates_dv->v + dsgemodel_ps->nlags;


   lht_dv = CreateVector_lf(dsgemodel_ps->fss);
   lhtst_dm = CreateMatrix_lf(smodel_ps->sv->nstates, dsgemodel_ps->fss);


   for (ti=dsgemodel_ps->nlags+1; ti<=dsgemodel_ps->nSample; ti++) //Note <= nSample because DW's function takes time as base-1.
   {
      tbase0i = ti - dsgemodel_ps->nlags - 1;

      //--- Getting logLH(y_t | Y_{t-1}, s_t, theta, q).
      for (si=smodel_ps->sv->nstates-1; si>=0; si--)
         lhtst_dm->M[mos(si,tbase0i,lhtst_dm->nrows)] = LogConditionalLikelihood(si,ti,smodel_ps);

      //--- Getting likelihood function at time t.
      lht_dv->v[tbase0i] = LogConditionalLikelihood_StatesIntegratedOut(ti, smodel_ps); //DW's function.
   }
   lhtst_dm->flag = M_GE;
   lht_dv->flag = V_DEF;



   //--- 1. Dates 1-by-fss.
   WriteVector(fptr_matlab2, &dates_effective_sdv, " %.16e ");
   //--- 2. Likelihood at time t 1-by-fss.
   fprintf(fptr_matlab2,"\n");
   WriteVector(fptr_matlab2, lht_dv, " %.16e ");
   //--- 3. Conditional (on s_t) likelihood at time 1 h-by-fss.
   fprintf(fptr_matlab2,"\n");
   WriteMatrix(fptr_matlab2, lhtst_dm, " %.16e ");
   //--- 4.ln(P(y[t] | Y[t-1], Z[t], s[t], theta, q)).
   fprintf(fptr_matlab2,"\n-------- ln(P(S[T] | theta, q)): -----------\n");
   fprintf(fptr_matlab2, " %.16e ", LogConditionalPrior_S(smodel_ps));


   //===
   DestroyVector_lf(lht_dv);
   DestroyMatrix_lf(lhtst_dm);
}
//---
static void ftd_FprintMatlab1(FILE *fptr_matlab1, struct TStateModel_tag *smodel_ps)
{
   int ti, tbase0i, si;
   //--- Taken out of Waggoner's model.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps = dsgemodel_ps->kalfilmsinputs_1stapp_ps;
   //--- Accessible variables.
   TSdvector dates_effective_sdv;  //fss-by-1 dates
   TSdcell *etdata_dc = kalfilmsinputs_1stapp_ps->etdata_dc; //ny-by-nst-by-T (T=fss+nlags), forecast errors e_t in the likelihood.
   TSdvector e_si_sdv; //ny-by-1.
   //===
   TSdvector et_sdv; //ny-by-1.
   TSdvector *lht_dv = NULL; //1-by-fss.
   TSdmatrix *et_dm = NULL; //ny-by-fss.

   dates_effective_sdv.flag = V_DEF;
   dates_effective_sdv.n = dsgemodel_ps->fss;
   dates_effective_sdv.v = dsgemodel_ps->dates_dv->v + dsgemodel_ps->nlags;


   lht_dv = CreateVector_lf(dsgemodel_ps->fss);
   et_dm = CreateMatrix_lf(etdata_dc->C[0]->nrows, dsgemodel_ps->fss);

   e_si_sdv.n = et_sdv.n = etdata_dc->C[0]->nrows;
   e_si_sdv.flag = et_sdv.flag = V_DEF;


   for (ti=dsgemodel_ps->nlags+1; ti<=dsgemodel_ps->nSample; ti++) //Note <= nSample because DW's function takes time as base-1.
   {
      tbase0i = ti - dsgemodel_ps->nlags - 1;

      //--- Getting forecasting errors at time t.
      et_sdv.v = et_dm->M + tbase0i*et_sdv.n;
      InitializeConstantVector_lf(&et_sdv, 0.0);
      for (si=smodel_ps->sv->nstates-1; si>=0; si--)
      {
         e_si_sdv.v = etdata_dc->C[ti-1]->M + si*e_si_sdv.n; //ti-1: base-0 counting for my own structure etdata_dc.
         ScalarTimesVectorUpdate(&et_sdv, ProbabilityStateConditionalPrevious(si, ti, smodel_ps), &e_si_sdv); //Prob*() is DW's function.
      }

      //--- Getting likelihood function at time t.
      lht_dv->v[tbase0i] = LogConditionalLikelihood_StatesIntegratedOut(ti, smodel_ps); //DW's function.
   }
   et_dm->flag = M_GE;
   lht_dv->flag = V_DEF;

   //--- 1. Dates 1-by-fss.
   WriteVector(fptr_matlab1, &dates_effective_sdv, " %.16e ");
   //--- 2. Likelihood at time t 1-by-fss.
   fprintf(fptr_matlab1,"\n");
   WriteVector(fptr_matlab1, lht_dv, " %.16e ");
   //--- 3. One-step forecast errors at time t ny-by-fss.
   fprintf(fptr_matlab1,"\n");
   WriteMatrix(fptr_matlab1, et_dm, " %.16e ");

   //===
   DestroyVector_lf(lht_dv);
   DestroyMatrix_lf(et_dm);
}
//---
static void ftd_FprintMatlab(FILE *fptr_matlab, struct TStateModel_tag *smodel_ps)
{
   //--- Taken out of Waggoner's model.
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   TSdvector dates_effective_sdv;  //fss-by-1 dates

   dates_effective_sdv.flag = V_DEF;
   dates_effective_sdv.n = dsgemodel_ps->fss;
   dates_effective_sdv.v = dsgemodel_ps->dates_dv->v + dsgemodel_ps->nlags;

   //--- 1. Dates 1-by-(fss+1).
   fprintf(fptr_matlab," %.16e ", dates_effective_sdv.v[0]-0.25); //Making it +1 from fss.
   WriteVector(fptr_matlab, &dates_effective_sdv, " %.16e ");
   //--- 2. Marginal posterior probabilities of base regimes. *-by-(fss+1).
   if (!dsgemodel_ps->indx_const)  //Regime-switching model.
   {
      fprintf(fptr_matlab,"\n");
      WriteCellTranspose(fptr_matlab, dsgemodel_ps->BaseProbS_dc, " %.16e ");
   }
}
//---
static void ftd_FprintEstResults(FILE *fptr_kernel7x_dv_localmodes, TSminpack *minpack_ps)
{
   TSdvector *x_dv = minpack_ps->x_dv;
   struct TStateModel_tag *smodel_ps = ((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->smodel_ps;
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;

   fprintf(fptr_kernel7x_dv_localmodes, " %0.16e ", dsgemodel_ps->peaklogpost);
   WriteVector(fptr_kernel7x_dv_localmodes, x_dv, " %0.16e ");
   fflush(fptr_kernel7x_dv_localmodes);
}
//---
static void ftd_FprintOutput(FILE *fptr_output, TSminpack *minpack_ps)
{
   static int header1=0;
   int si;
   //---
   TSdvector *x_dv = minpack_ps->x_dv;
   TSdvector *g_dv = minpack_ps->g_dv;
   struct TSetc_csminwel_tag *etc_csminwel_ps = (TSetc_csminwel *)minpack_ps->etc_package_ps;
   struct TSargs_blockcsminwel_tag *args_blockcsminwel_ps = ((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->args_blockcsminwel_ps;
   struct TStateModel_tag *smodel_ps = ((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->smodel_ps;
   struct TSdsgemodel_tag *dsgemodel_ps = (struct TSdsgemodel_tag *)smodel_ps->theta;
   struct TSfinalpars_tag *finalpars_ps = dsgemodel_ps->finalpars_ps;
   struct TSgensys_tag *gensys_ps = dsgemodel_ps->gensys_ps;
   struct TSDWRegimeControlVariables_tag *DWRegimeControlVariables_ps = dsgemodel_ps->DWRegimeControlVariables_ps;
   int *n_st_var = DWRegimeControlVariables_ps->n_st_var;
   int *n_st_coef = DWRegimeControlVariables_ps->n_st_coef;
   //
   int fss = dsgemodel_ps->fss;
   //--- Pointers only.
   TSdvector modelx_sdv;
   TSdvector modelx0_sdv;
   TSdvector transition_qs_sdv;
   TSdvector transition_q0s_sdv;

   modelx_sdv.n = modelx0_sdv.n = dsgemodel_ps->xphi_dv->n;
   modelx_sdv.v = x_dv->v;
   modelx0_sdv.v = minpack_ps->x0_dv->v;

   dsgemodel_ps->program_hours = difftime(dsgemodel_ps->prog_endtime,dsgemodel_ps->prog_begtime)/60.0/60.0;
   
   if (header1==0) {
      fprintf(fptr_output, "===================== Reporting final results =====================\n");
      header1=1;        //After this reset, when this function is called second time, this header will not be printed again.
   }

   
   
   //--- Printing out model TAG.
   fprintf(fptr_output, "\n----------------- Model tag is: %s -----------------\n", MODELTAG);
   fprintf(fptr_output, "\n----------------- Model restriction type: indx_hard_restrict = %d -----------------\n", dsgemodel_ps->indx_hard_restrict);
   fprintf(fptr_output, "\n----------------- Number of lags: nlags = %d -----------------\n", dsgemodel_ps->nlags);
   fprintf(fptr_output, "\n----------------- Effective sample size (excluding lags): fss = %d -----------------\n", dsgemodel_ps->fss);
   fprintf(fptr_output, "\n----------------- Sample period (including lags) is -----------------\n");
   WriteVector_int(fptr_output, dsgemodel_ps->BeginEndSample_iv);
   
   
   //--------------------------- Printing some key flags. ---------------------------
   fprintf(fptr_output, "indxEstFinder = %d -- if 1, estimation takes place; otherwise, no estimation.\n", dsgemodel_ps->indxEstFinder);
   fprintf(fptr_output, "\n\n");


   //--------------------------- Printing final outputs. ---------------------------
   fprintf(fptr_output, "Log of the scale for the prior density that must satisfy certain restrictions: %.5e.\n\n", dsgemodel_ps->scale4logpriordensity);
   fprintf(fptr_output, "Gensys existence and uniqueness at the posterior estimate: [%d, %d].\n", gensys_ps->eu_iv->v[0], gensys_ps->eu_iv->v[1]);
   fprintf(fptr_output, "Starting method:  0 -- continuing from the last estimated results contained in filename_sp;\n"
                        "                  1 -- starts from the fixed values for xphi_dv, manually keyed in datainpu_setup.prn;\n"
                        "                  2 -- randomly or arbitarily selects the initial starting values for the MLE or posterior estimate.\n");
   fprintf(fptr_output, "Starting value is %d.  If the method is indexed by 2, the random seed number is %d\n", dsgemodel_ps->indxStartValuesForMin, dsgemodel_ps->randomseed);
   //fprintf(fptr_output, "Index for the MLE is %d (0: posterior estimates; 1: MLE with no prior)\n", constmodpars_ps->indxMLE);
   fprintf(fptr_output, "Total computing hours for the whole program: %10.3f\n", dsgemodel_ps->program_hours);

   //=== Minimization problem.
   fprintf(fptr_output, "\n------- Output arguments from the csminwel package. -------\n");
   fprintf(fptr_output, "Numerical gradient associated with the posterior estimates g_dv:\n");
   WriteVector(fptr_output, g_dv, AFORMAT);
   if (!dsgemodel_ps->indx_const)
   {
      fprintf(fptr_output, "The end convergence criterion for the function value: %.5e\n", args_blockcsminwel_ps->criterion_end);
      fprintf(fptr_output, "Initial value for the diagonal of inverse Hessian in the quasi-Newton search: %.5e\n", args_blockcsminwel_ps->ini_h_scale);
      fprintf(fptr_output, "Blockwise step sizes of non-automatic portion of numerical gradient:\n");
      WriteVector(fptr_output, args_blockcsminwel_ps->gradstps_csminwel_dv, " %.5e ");
      //+
      transition_qs_sdv.flag = transition_q0s_sdv.flag = V_DEF;
      transition_qs_sdv.n = transition_q0s_sdv.n = NumberFreeParametersQ(smodel_ps);
      transition_qs_sdv.v = x_dv->v + modelx_sdv.n;
      transition_q0s_sdv.v = minpack_ps->x0_dv->v + modelx0_sdv.n;
   }
   else  //Constant-parameter case
   {
      fprintf(fptr_output, "Overall convergence criterion for the function value: %.5e\n", etc_csminwel_ps->crit);
      fprintf(fptr_output, "Initial value for the diagonal of inverse Hessian in the quasi-Newton search: %.5e\n", etc_csminwel_ps->ini_h_csminwel);
      fprintf(fptr_output, "Step size of non-automatic portion of numerical gradient: %.5e\n", etc_csminwel_ps->gradstps_csminwel);
   }
   fprintf(fptr_output, "Method of numerical gradient is: %d, where 1 represents forward difference and 2 represents central difference\n", etc_csminwel_ps->indxnumgrad_csminwel);
   fprintf(fptr_output, "Number of iterations taken by csminwel: %d\n", etc_csminwel_ps->niter);
   fprintf(fptr_output, "Number of function evaluations used by csminwel: %d\n", etc_csminwel_ps->fcount);
   fprintf(fptr_output, "Return code for the terminating condition: %d\n", etc_csminwel_ps->retcode);
   //+
   fprintf(fptr_output, "\n\n--------------\n");
   fprintf(fptr_output, "Initial or starting values of free model parameters x0_dv (excluding Q elements):\n");
   WriteVector(fptr_output, &modelx0_sdv, AFORMAT);
   if (!dsgemodel_ps->indx_const)  //Markov-switching case.
   {
      fprintf(fptr_output, "Initial or starting values of free Q elements transition_q0s_sdv.v:\n");
      WriteVector(fptr_output, &transition_q0s_sdv, AFORMAT);
   }
   fprintf(fptr_output, "The log value of posterior kernel at the initial values: %.16e\n\n", -minpack_ps->fret0);
   //+
   fprintf(fptr_output, "Posterior estimates of free model parameters x_dv (excluding Q elements):\n");
   WriteVector(fptr_output, &modelx_sdv, AFORMAT);
   if (!dsgemodel_ps->indx_const)  //Markov-switching case.
   {
      fprintf(fptr_output, "Posterior estimates of free Q elements transition_qs_sdv.v:\n");
      WriteVector(fptr_output, &transition_qs_sdv, AFORMAT);
   }
   fprintf(fptr_output, "------- Note that the log posterior kernel below must be adjusted by adding the log scale factor and that the log scale factor may.\n");
   fprintf(fptr_output, "-------   vary for a given set of parameters because of MCMC rounding errors in calculating scale4logpriordensity for the prior.\n");
   fprintf(fptr_output, "Peak value of logPosterior: %.16e\n", dsgemodel_ps->peaklogpost);
   fprintf(fptr_output, "log scale factor from the prior to be added to log posterior kernel: %.16e\n", dsgemodel_ps->scale4logpriordensity);   
   fprintf(fptr_output, "LH (pdf) value at the posterior peak: %.16e\n", dsgemodel_ps->LHvalatpostpeak);


   //=== Final individual parameters.
   #include CFprintOutput_TAG


   fprintf(fptr_output, "\n%%--- Transition probibility matrices for different regime variables: ---\n");
   for (si=0; si<smodel_ps->sv->n_state_variables; si++)
   {
      fprintf(fptr_output, "Q(:,:,%d) = [\n", si+1);
      WriteMatrix(fptr_output, smodel_ps->sv->state_variable[si]->baseQ, " %.7g ");
      fprintf(fptr_output, "];\n\n");
   }


   //--- Printing out covariance matrix at the log posterior peak for Metropolis use or asymptotic reports.
   if (dsgemodel_ps->Omega_dm->flag)
   {
      fprintf(fptr_output, "\n%%--- Covariance matrix at the log posterior peak ---\n");
      fprintf(fptr_output, "CovMatrix = [\n");
      WriteMatrix(fptr_output, dsgemodel_ps->Omega_dm, " %.12e ");
      fprintf(fptr_output, "];\n\n");
   }
   //--- Printing out Hessian at the log posterior peak and its inverse (Omega) is used for Metropolis use or asymptotic reports.
   if (dsgemodel_ps->Hessian_dm->flag)
   {
      fprintf(fptr_output, "\n%%--- Hessian matrix at the log posterior peak ---\n");
      fprintf(fptr_output, "Hessian = [\n");
      WriteMatrix(fptr_output, dsgemodel_ps->Hessian_dm, " %.12e ");
      fprintf(fptr_output, "];\n\n");
   }
}
#undef RFORMAT
#undef NFORMAT
#undef AFORMAT
#undef ARFORMAT
