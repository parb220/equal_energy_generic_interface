#ifndef _MPI_PARAMETER_HEADER_
#define _MPI_PARAMETER_HEADER_

const size_t N_MESSAGE = 1000; // to be safe
const int LENGTH_INDEX = 0; 
const int LEVEL_INDEX = 1; 
const int BURN_INDEX = 2; 
const int thin_INDEX = 3;
const int THIN_INDEX = 4; 
const int GROUP_INDEX = 5; 
const int PEE_INDEX = 6; 
const int RETURN_INDEX_1 = 7; 
const int RETURN_INDEX_2 = 8; 
const int RETURN_INDEX_3 = 9; 
const int RETURN_INDEX_4 = 10; 
const int RETURN_INDEX_5 = 11; 
const int SCALE_INDEX = 12;
const int KPLUS_INDEX = 14;
const int K_INDEX = 14;
const int RESERVE_INDEX = 100; 

const int BURN_IN_LENGTH = 5000; 
const int SIMULATION_LENGTH = 50000; 

const int TUNE_TAG = 19; 
const int TUNE_TAG_DIAG = 20;
const int TUNE_TAG_VARIANCE = 21;  
const int TUNE_TAG_VARIANCE_BASED_ON_DIAG = 25; 
const int TUNE_TAG_VARIANCE_BASED_ON_HIGHER = 23;
const int TUNE_TAG_VARIANCE_BASED_ON_LOWER = 26; 
const int TUNE_TAG_VARIANCE_BASED_ON_CURRENT = 24; 
const int TUNE_TAG_BEFORE_SIMULATION = 28; 
const int TUNE_TAG_AFTER_SIMULATION = 29; 

const int TUNE_TAG_DIAG_RANDOM_BLOCK = 80;
const int TUNE_TAG_VARIANCE_RANDOM_BLOCK = 81;  
const int TUNE_TAG_VARIANCE_BASED_ON_DIAG_RANDOM_BLOCK = 85; 
const int TUNE_TAG_VARIANCE_BASED_ON_HIGHER_RANDOM_BLOCK = 83;
const int TUNE_TAG_VARIANCE_BASED_ON_LOWER_RANDOM_BLOCK = 86; 
const int TUNE_TAG_VARIANCE_BASED_ON_CURRENT_RANDOM_BLOCK = 84; 

const int TUNE_TAG_SIMULATION_FIRST = 32; 
const int TUNE_TAG_SIMULATION_FIRST_MULTIPLE_BLOCK = 37; 
const int TUNE_TAG_SIMULATION_SECOND = 42;
const int TUNE_TAG_SIMULATION_SECOND_MULTIPLE_BLOCK = 47;
const int TUNE_TAG_INITIAL = 48;                             // added 9/19/2014

const int TRACKING_TAG_SIMULATION = 50;
const int TRACKING_TAG_SIMULATION_MULTIPLE_BLOCK = 51; 

const int SIMULATION_TAG = 60; 
const int NP_SOL_TAG = 61; 
const int SIMULATION_TAG_MULTIPLE_BLOCK = 62; 
const int SIMULATION_TAG_INITIAL_METROPOLIS = 63;                       // added 9/19/2014
const int SIMULATION_TAG_INITIAL_INDEPENDENT = 64;          // added 9/26/2014

const int HILL_CLIMB_TAG = 70; 
const int GMM_SIMULATION_TAG = 71; 

const int BINNING_INFO  = 1000; 
 
const int END_TAG = 0;  
#endif
