#ifndef _MPI_PARAMETER_HEADER_
#define _MPI_PARAMETER_HEADER_

const size_t N_MESSAGE = 1000; // to be safe
const int LENGTH_INDEX = 0; 
const int LEVEL_INDEX = 1; 
const int BURN_INDEX = 2; 
const int SCALE_INDEX = 3; 
const int THIN_INDEX = 4; 
const int GROUP_INDEX = 5; 
const int PEE_INDEX = 6; 
const int RETURN_INDEX_1 = 7; 
const int RETURN_INDEX_2 = 8; 
const int RETURN_INDEX_3 = 9; 
const int RETURN_INDEX_4 = 10; 
const int RETURN_INDEX_5 = 11; 
const int RESERVE_INDEX = 100; 
const int GROUP_NUMBER_INDEX = 101; 

const int BURN_IN_LENGTH = 5000; 
const int SIMULATION_LENGTH = 50000; 

const int TUNE_TAG_BEFORE_SIMULATION = 28; 

const int TUNE_TAG_SIMULATION_FIRST = 32; 

const int SIMULATION_TAG = 60; 
const int SIMULATION_PRIOR_TAG = 63; 

const int BINNING_INFO  = 1000; 

const int SCALE_MATRIX_FIT_TAG = 91; 

const int END_TAG = 0;  
const int TASK_LENGTH_LONG = 1000; 
const int TASK_LENGTH_SHORT = 100; 
#endif
