#ifndef __DSGELINV1_CEXPS_H__

   #define __DSGELINV1_CEXPS_H__

   //#define _IMSLOPTPACKAGE_


   //============= Types of hard restrictions on free parameters: ====================
   //======= Benchmark model (TAG): mod_6v. =======
   //=== mod: 6 variables, no decpreciation shock, gtheta derived, and by fixed.
   //---   mod0: default model (persistence for all exogenous shocks)
   //---   mod1: restricting glambdastar and geta but everything else is the same as model 0.
   //---   mod2: restricting gbeta but everything else is the same as model 0.
   //---   mod3: restricting galpha and geta with beta derived, but everything else is the same as model 0.  USED as checking modc203
   //---   mod4: restricting galpha, glambdastar, and geta but everything else is the same as model 0.
   //---   mod5: restricting galpha and geta with beta derived -- the same as model 3 but with R to be freely estimated.
   //---   mod6: restricting galpha and geta, -- the same as model 3 but with R to be freely estimated.
   //=== modb: 6 variables, depreciation shock, gtheta fixed, and by derived.
   //---   modb103: adding depreciation shock and deleting the debt series (5 variables total),
   //                  but everything else is the same as model mod3 (which has the debt series).  OK.
   //---   modb105: adding depreciation shock and deleting the debt series (5 variables total),
   //                  but everything else is the same as model mod5 (which has the debt series).  OK.
   //=== modc: 6 variables, no depreciation shock, gtheta fixed, and by derived.
   //---   modc203: no depreciation shock and exactly the same as mod3.  USED FOR THE PAPER.
   //---   modc205: no depreciation shock and exactly the same as mod5.  USED as checking.
   //======= Alternative model (TAG): modalt, suggested by referee, where the credit constraint depends only on capital value. =======   
   //=== modc_alt: 
   
   #define __BENCHMARK__
   #undef __CAPITALCONSTRAINTONLY__

   #if defined(__BENCHMARK__)
   
      #include "dsgelinv1_benchmark.h"
      //--- Automatically exported from ConvertMatlab2C.py and from Mathematica.
      #define CInputGensysForm_TAG  "Cinput_Gensys_All.cexps"            //Called by dsgelinv1_comfuns.c
      #define CInputMeasureMatrices_TAG  "Cinput_MeasureEqns_All.cexps"  //Called by dsgelinv1_comfuns.c
      #define CInputDerivedSSPars_TAG  "CInput_DerivedPars_ALL.cexps"   //Called by dsgelinv1_MODEL.c, where MODEL = "benchmark," for example
      //--- Manually keyed in.
      #define CRefresh_Kalman_allcases_TAG  "CRefresh_Kalman_All.cexps"         //Called by dsgelinv1_comfuns.c
      #define CFprintOutput_TAG  "CFprintOutput_ALL.cexps"                          //Called by dsgelinv1_estmcmc.c
   #elif defined(__MODEL2__)

      #include "dsgelinv1_mod2_3v_const.h"
      //--- Automatically exported from ConvertMatlab2C.py and from Mathematica.
      #define CInputGensysForm_TAG  "CInputGensysForm_mod1_3v_const.cexps" //mod2 and mod1 same.
      #define CInputMeasureMatrices_TAG  "CInputMeasureMatrices_mod1_3v_const.cexps" //mod2 and mod1 same.
      #define CInputDerivedSSPars_TAG  "CInputDerivedSSPars_mod2_3v_const.cexps"
      //--- Manually keyed in.
      #define CRefresh_kalfilms_allcases_TAG  "CRefresh_kalfilms_allcases_mod1_3v_const.cexps"
                        //mod2 and mod1 same.
      #define CFprintOutput_TAG  "CFprintOutput_mod2_3v_const.cexps"
   #endif



#endif
