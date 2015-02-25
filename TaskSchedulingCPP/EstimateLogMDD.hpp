#ifndef _ESTIMATE_LOG_MDD_HEADER_
#define _ESTIMATE_LOG_MDD_HEADER_

const int POSTERIOR_HEATED = 0;
const int LIKELIHOOD_HEATED = 1;
const int PRIOR_ONLY = 2;


bool compare_CSampleIDWeight_BasedOnEnergy(const CSampleIDWeight &i, const CSampleIDWeight &j);
bool compare_CSampleIDWeight_BasedOnID(const CSampleIDWeight &i, const CSampleIDWeight &j);

// LogMDD using Gaussian as proposal
double EstimateLogMDD_gaussian(CEquiEnergyModel &model, int stage, int posterior_type=POSTERIOR_HEATED);
double LogMDD_gaussian(const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, double t, int posterior_type = POSTERIOR_HEATED);
vector<double>EstimateLogMDD_gaussian(CEquiEnergyModel &model, int stage, int nGroup, int posterior_type=POSTERIOR_HEATED);

// LogMDD using elliptical as proposal
double EstimateLogMDD(CEquiEnergyModel &model, int stage, int proposal_type, int posterior_type=POSTERIOR_HEATED);
double LogMDD(const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, double t, int proposal_type, int posterior_type = POSTERIOR_HEATED);
vector<double> WeightFromGaussianSample(int N, const TDenseVector &center, const TDenseMatrix &scale, CEquiEnergyModel &model, double temperature, int posterior_type=POSTERIOR_HEATED);
vector<double>EstimateLogMDD(CEquiEnergyModel &model, int stage, int proposal_type, int nGroup, int posterior_type=POSTERIOR_HEATED);
double CheckLogMDDConvergency(vector<CSampleIDWeight> &sample, CEquiEnergyModel &model, int stage, double t, int proposal_type, double &average_logMDD, double &std_logMDD, int posterior_type = POSTERIOR_HEATED);

// LogMDD using samples from another stage as proposal
double EstimateLogMDD(CEquiEnergyModel &model, int stage, int previous_stage, double logMDD_previous, int posterior_type=POSTERIOR_HEATED);
double LogMDD(const vector<CSampleIDWeight> &proposal, const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, int previous_stage, int current_stage, double logMDD_previous, int posterior_type=POSTERIOR_HEATED);

// LogMDD with importance sampling
double LogMDD_Importance(const vector<CSampleIDWeight> &proposal, CEquiEnergyModel &model,  int previous_stage, int current_stage, double logMDD_previous, int posterior_type=POSTERIOR_HEATED);

#endif 
