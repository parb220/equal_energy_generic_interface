#include <sstream>
#include <string>
#include "CEquiEnergyModel.h"
#include "CSampleIDWeight.h"
#include "CStorageHead.h"
#include "dw_matrix.h"
#include "storage_parameter.h" 
#include "prcsn.h"
#include "dw_math.h"
#include "mdd_function.h"

using namespace std; 

void EstimateLogMDD(CEquiEnergyModel &model, int level, double logMDD_previous_bridge, double logMDD_previous_mueller, double &logMDD_bridge, double &logMDD_mueller)
{
	vector<CSampleIDWeight> proposal, posterior; 
	if (!model.storage->DrawAllSample(level+1, proposal) || !model.storage->DrawAllSample(level, posterior)) 
	{
		cerr << "EstimateLogMDD:: error occurred when loading all samples.\n"; 
		abort(); 
	}

	// dimension 0: proposal density
	// dimension 1: posterior kernel
	TMatrix proposal_value_bridge=CreateMatrix(proposal.size(), 2); 
	TMatrix posterior_value_bridge=CreateMatrix(posterior.size(), 2);  
	TMatrix proposal_value_mueller=CreateMatrix(proposal.size(), 2); 
	TMatrix posterior_value_mueller=CreateMatrix(posterior.size(), 2); 

	if (level == (int) model.parameter->number_energy_level-1)
	{
		stringstream convert;
                convert.str(string());
                convert <<  model.parameter->run_id << "/" << model.parameter->run_id << GM_MEAN_COVARIANCE;
               	string gm_file = model.parameter->storage_dir + convert.str();
                if (!model.ReadGaussianMixtureModelParameters(gm_file) )
                {
                	cerr << "Error occurred while reading Gaussian mixture model parameters from " << gm_file << endl;
                	abort();
                }
		for (int i=0; i<(int)proposal.size(); i++)
			ElementM(proposal_value_bridge, i, 0) = ElementM(proposal_value_mueller, i, 0) =  model.GMM_LogPDF(proposal[i]); 
		for (int i=1; i<(int)posterior.size(); i++)
			ElementM(posterior_value_bridge, i, 0) = ElementM(posterior_value_mueller, i, 0) =  model.GMM_LogPDF(posterior[i]);  
	}
	else 
	{
		for(int i=0; i<(int)proposal.size(); i++)
		{
			ElementM(proposal_value_bridge, i, 0) = proposal[i].weight/model.parameter->t[level+1]-logMDD_previous_bridge; 
			ElementM(proposal_value_mueller, i, 0) = proposal[i].weight/model.parameter->t[level+1]-logMDD_previous_mueller; 
		}
		for(int i=0; i<(int)posterior.size(); i++)
		{
			ElementM(posterior_value_bridge, i, 0) = posterior[i].weight/model.parameter->t[level+1]-logMDD_previous_bridge; 
			ElementM(proposal_value_mueller, i, 0) = posterior[i].weight/model.parameter->t[level+1]-logMDD_previous_mueller; 
		}	
	}
	
	for(int i=0; i<(int)proposal.size(); i++)
        	ElementM(proposal_value_bridge, i, 1) = ElementM(proposal_value_mueller, i, 1) = proposal[i].weight/model.parameter->t[level];
        for(int i=0; i<(int)posterior.size(); i++)
                ElementM(posterior_value_bridge, i, 1) = ElementM(posterior_value_mueller, i, 1) = posterior[i].weight/model.parameter->t[level];

	logMDD_bridge = ComputeLogMarginalDensity_Bridge(proposal_value_bridge, posterior_value_bridge); 
	int in_P1, in_P2; 
	logMDD_mueller = ComputeLogMarginalDensity_Mueller(proposal_value_mueller, posterior_value_mueller, &in_P1, &in_P2); 
	
}

/*
double ComputeLogMDD_Bridge(const TDenseMatrix &proposal, const TDenseMatrix &posterior)
{
	double min_c, max_c, mid_c=0.0, diff;

	// Bracket the zero
	if ((diff=BridgeDifference(proposal,posterior,mid_c)) < 0.0)
    	{
      	max_c=mid_c;
      	for (min_c=-1.0; min_c > MINUS_INFINITY; max_c=min_c, min_c*=10)
        	if ((diff=BridgeDifference(proposal,posterior,min_c)) > 0)
			break;
      	if (min_c <= MINUS_INFINITY) 
		return min_c;
    	}
  	else
    	{
      		min_c=mid_c;
      		for (max_c=1.0; max_c < PLUS_INFINITY; min_c=max_c, max_c*=10)
        		if ((diff=BridgeDifference(proposal,posterior,max_c)) < 0)
				break;
      		if (max_c >= PLUS_INFINITY) 
			return max_c;
    	}

	// Divide and conququer
	diff=BridgeDifference(proposal,posterior,mid_c=(min_c + max_c)/2.0);
  	for (int i=0; i < 300; i++)
    	{
      		if (diff > 0)
        		min_c=mid_c;
      		else
        		max_c=mid_c;
      		if ((fabs(diff=BridgeDifference(proposal,posterior,mid_c=(min_c + max_c)/2.0)) < SQRT_MACHINE_EPSILON)) 
			break;
    	}
  	return mid_c;
}

double BridgeDifference(const TDenseMatrix &proposal, const TDenseMatrix &posterior, double logr)
{	
	double x, sum1=MINUS_INFINITY, sum2=MINUS_INFINITY; 
 	int n1=posterior.rows, n2=proposal.rows;

	double r = (double)n2/(double)n1; 
  	for (int i=n2-1; i >= 0; i--)
    		if ((x=proposal(i,0)+logr-proposal(i,1)) < 0)
      			sum2=AddLogs(sum2,-log(1+r*exp(x)));
    		else
      			sum2=AddLogs(sum2,-x-log(exp(-x)+r));

	r=(double)n1/(double)n2; 
  	for (int i=n1-1; i >= 0; i--)
    		if ((x=posterior(i,1)-posterior(i,0)-logr) < 0)
      			sum1=AddLogs(sum1,-log(1+r*exp(x)));
    		else
      			sum1=AddLogs(sum1,-x-log(exp(-x)+r));

  	return sum2 - sum1;
}
*/


