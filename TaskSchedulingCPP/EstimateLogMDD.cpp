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
#include "mdd.hpp"
#include "dw_elliptical.h"

using namespace std; 

double EstimateLogMDD(CEquiEnergyModel &model, int level, int previous_level, double logMDD_previous); 
double EstimateLogMDD(CEquiEnergyModel &model, int level, int proposal_type);
double EstimateLogMDD_gaussian(CEquiEnergyModel &model, int level); 

TMatrix CreateProposalMatrix(int ndraws, CEquiEnergyModel &model, const TElliptical *elliptical); 
TMatrix CreatePosteriorMatrix(const vector<CSampleIDWeight> &sample, const TElliptical *elliptical);
TMatrix CreateGaussianProposalMatrix(int ndraws, CEquiEnergyModel &model, const TDenseVector &center, const TDenseMatrix &scale);
TMatrix CreateGaussianPosteriorMatrix(const vector<CSampleIDWeight> &sample, const TDenseVector &center, const TDenseMatrix &scale);
TMatrix CreateProposalMatrix(int ndraws, CEquiEnergyModel &model, double t, const TElliptical *elliptical); 
TMatrix CreatePosteriorMatrix(const vector<CSampleIDWeight> &sample, double t, const TElliptical *elliptical);
TMatrix CreateGaussianProposalMatrix(int ndraws, CEquiEnergyModel &model, double t, const TDenseVector &center, const TDenseMatrix &scale);
TMatrix CreateGaussianPosteriorMatrix(const vector<CSampleIDWeight> &sample, double t, const TDenseVector &center, const TDenseMatrix &scale);
TDenseVector GetRadiusFromSample(const vector<CSampleIDWeight> &sample, const TDenseVector &center, const TDenseMatrix &scale); 
bool GetCenterScaleFromSample(const vector<CSampleIDWeight> &sample, TDenseVector &center, TDenseMatrix &scale); 

TMatrix CreateGaussianPosteriorMatrix(const vector<CSampleIDWeight> &sample, const TDenseVector &center, const TDenseMatrix &scale)
{
	TMatrix posterior = CreateMatrix(sample.size(), 2); 
	double logdensity; 
	TDenseVector y(sample[0].data.dim, 0.0); 
	TDenseMatrix scale_inverse = Inverse(scale); 
	
	for (int i=0; i<(int)sample.size(); i++)
	{
		y = scale_inverse*(sample[i].data-center);  // y = Inverse(scale)*(x-center); 
		logdensity = -0.918938533204673 -0.5*y[0]*y[0];
		for (int j=1; j<y.dim; j++)
                        logdensity += -0.918938533204673 -0.5*y[j]*y[j];
		ElementM(posterior,i,0)=logdensity;
		ElementM(posterior, i, 1) = sample[i].weight; 
	}
	return posterior; 
}

TMatrix CreateGaussianPosteriorMatrix(const vector<CSampleIDWeight> &sample, double t, const TDenseVector &center, const TDenseMatrix &scale)
{
	TMatrix posterior = CreateMatrix(sample.size(), 2); 
	double logdensity; 
	TDenseVector y(sample[0].data.dim, 0.0); 
	TDenseMatrix scale_inverse = Inverse(scale); 
	
	for (int i=0; i<(int)sample.size(); i++)
	{
		y = scale_inverse*(sample[i].data-center);  // y = Inverse(scale)*(x-center); 
		logdensity = -0.918938533204673 -0.5*y[0]*y[0];
		for (int j=1; j<y.dim; j++)
                        logdensity += -0.918938533204673 -0.5*y[j]*y[j];
		ElementM(posterior,i,0)=logdensity;
		ElementM(posterior, i, 1) = sample[i].weight/t; 
	}
	return posterior; 
}

TMatrix CreateGaussianProposalMatrix(int ndraws, CEquiEnergyModel &model, const TDenseVector &center, const TDenseMatrix &scale)
{
	if (ndraws <= 0)
		return (TMatrix)NULL; 
	
	TMatrix proposal = CreateMatrix(ndraws, 2); 
	TDenseVector x(center.dim,0.0), y(center.dim, 0.0); 
	double log_density, log_posterior; 
	for (int i=0; i<ndraws; i++)
	{
		x.RandomNormal(); 
		log_density = -0.918938533204673 -0.5*x[0]*x[0];
		for (int j=1; j<x.dim; j++)
			log_density += -0.918938533204673 -0.5*x[j]*x[j];

		y = scale * x + center; // y = center + scale *x; 
		log_posterior = model.log_posterior_function(y.vector, y.dim); 
		ElementM(proposal, i, 0) = log_density;
                ElementM(proposal, i, 1) = log_posterior;
	}
	return proposal; 
}

TMatrix CreateGaussianProposalMatrix(int ndraws, CEquiEnergyModel &model, double t, const TDenseVector &center, const TDenseMatrix &scale)
{
	if (ndraws <= 0)
		return (TMatrix)NULL; 
	
	TMatrix proposal = CreateMatrix(ndraws, 2); 
	TDenseVector x(center.dim,0.0), y(center.dim, 0.0); 
	double log_density, log_posterior; 
	for (int i=0; i<ndraws; i++)
	{
		x.RandomNormal(); 
		log_density = -0.918938533204673 -0.5*x[0]*x[0];
		for (int j=1; j<x.dim; j++)
			log_density += -0.918938533204673 -0.5*x[j]*x[j];

		y = scale * x + center; // y = center + scale *x; 
		log_posterior = model.log_posterior_function(y.vector, y.dim); 
		ElementM(proposal, i, 0) = log_density;
                ElementM(proposal, i, 1) = log_posterior/t;
	}
	return proposal; 
}

TMatrix CreateProposalMatrix(int ndraws, CEquiEnergyModel &model, const TElliptical *elliptical)
{
	if (ndraws <= 0)
		return (TMatrix)NULL; 

	TMatrix proposal = CreateMatrix(ndraws, 2); 
	int dim = elliptical ? elliptical->dim : 0;
	double *draw = dim>0 ? new double[dim] : NULL; 

	double log_density, log_posterior;
	for (int i=0; i<ndraws; i++)
	{
		log_density = draw ?  LogDensityElliptical_Radius(DrawElliptical(draw, (TElliptical *)elliptical), (TElliptical *)elliptical) : 0.0;
		ElementM(proposal, i, 0) = log_density; 

		log_posterior = model.log_posterior_function(draw, dim); 
		ElementM(proposal, i, 1) = log_posterior; 
	}

	if (draw)
		delete [] draw; 
	return proposal; 
}

TMatrix CreateProposalMatrix(int ndraws, CEquiEnergyModel &model, double t, const TElliptical *elliptical)
{
	if (ndraws <= 0)
		return (TMatrix)NULL; 

	TMatrix proposal = CreateMatrix(ndraws, 2); 
	int dim = elliptical ? elliptical->dim : 0;
	double *draw = dim>0 ? new double[dim] : NULL; 

	double log_density, log_posterior;
	for (int i=0; i<ndraws; i++)
	{
		log_density = draw ?  LogDensityElliptical_Radius(DrawElliptical(draw, (TElliptical *)elliptical), (TElliptical *)elliptical) : 0.0;
		ElementM(proposal, i, 0) = log_density; 

		log_posterior = model.log_posterior_function(draw, dim); 
		ElementM(proposal, i, 1) = log_posterior/t; 
	}

	if (draw)
		delete [] draw; 
	return proposal; 
}

TMatrix CreatePosteriorMatrix(const vector<CSampleIDWeight> &sample, const TElliptical *elliptical)
{
	TMatrix posterior=CreateMatrix(sample.size(),2); 
	int dim = elliptical ? elliptical->dim:0; 

	double log_density; 
	for (int i=0; i<(int)sample.size(); i++)
	{
		log_density = dim> 0 ? LogDensityElliptical_Draw(sample[i].data.vector, (TElliptical *)elliptical) : 0.0; 
		ElementM(posterior, i, 0) = log_density; 
		ElementM(posterior, i, 1) = sample[i].weight; 
	}
	return posterior; 
}

TMatrix CreatePosteriorMatrix(const vector<CSampleIDWeight> &sample, double t, const TElliptical *elliptical)
{
	TMatrix posterior=CreateMatrix(sample.size(),2); 
	int dim = elliptical ? elliptical->dim:0; 

	double log_density; 
	for (int i=0; i<(int)sample.size(); i++)
	{
		log_density = dim> 0 ? LogDensityElliptical_Draw(sample[i].data.vector, (TElliptical *)elliptical) : 0.0; 
		ElementM(posterior, i, 0) = log_density; 
		ElementM(posterior, i, 1) = sample[i].weight/t; 
	}
	return posterior; 
}


TDenseVector GetRadiusFromSample(const vector<CSampleIDWeight> &sample, const TDenseVector &center, const TDenseMatrix &scale)
{
	TDenseVector Radii(sample.size(), 0.0); 
	TDenseMatrix quadratic_form = MultiplyTranspose(scale, scale); 
	quadratic_form.Inverse();
	quadratic_form = (quadratic_form + Transpose(quadratic_form)) * 0.5; // force symmetry
	
	TDenseVector theta; 
	for (int i=0; i<(int)sample.size(); i++)
	{
		theta = sample[i].data - center; 
		Radii[i] =  sqrt(InnerProduct(theta, theta, quadratic_form)); 
	}
	return Radii; 
}

bool GetCenterScaleFromSample(const vector<CSampleIDWeight> &sample, TDenseVector &center, TDenseMatrix &scale)
{
	if (sample.empty())
	{
		cerr << "GetCenterScaleFromSample() : sample empty." << endl; 
		return false; 
	}

	double max; 
	CSampleIDWeight mode; 
	TDenseVector sample_sum(sample[0].data.dim, 0.0); 
	TDenseMatrix sample_square(sample[0].data.dim, sample[0].data.dim, 0.0); 
	for (int i=0; i<(int)sample.size(); i++)
	{
		if (i==0 || sample[i].weight > max)
		{
			max = sample[i].weight; 
			mode = sample[i];  
		}
		sample_sum = sample_sum + sample[i].data; 
		sample_square = sample_square + OuterProduct(sample[i].data, sample[i].data); 
	}
	
	sample_sum = sample_sum * (1.0/(double)sample.size()); 
	sample_square = sample_square * (1.0/(double)sample.size()); 	

	center = mode.data; 
	sample_square = sample_square + OuterProduct(mode.data, mode.data); 
	sample_square = sample_square - OuterProduct(sample_sum, mode.data) - OuterProduct(mode.data, sample_sum); 
	sample_square = (sample_square+Transpose(sample_square))*0.5; 
	scale = Cholesky(sample_square); 
	scale = Transpose(scale); 
	return true; ; 
}

double EstimateLogMDD_gaussian(CEquiEnergyModel &model, int level)
{
	vector<CSampleIDWeight > posterior;
        bool unstructured = true;
        int data_size = model.current_sample.GetSize_Data();
        if (!model.storage->DrawAllSample(level, posterior, unstructured, data_size))
        {
                cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
                abort();
        }
	
	TDenseVector center(posterior[0].data.dim, 0.0);
        TDenseMatrix scale(posterior[0].data.dim, posterior[0].data.dim, 0.0);
        if(!GetCenterScaleFromSample(posterior, center, scale))
        {
                cerr << "Error occurred in GetCenterScaleFromSample()\n";
                abort();
        }

	TMatrix proposal_matrix = CreateGaussianProposalMatrix(posterior.size(), model, model.parameter->t[level], center, scale); 
	TMatrix posterior_matrix = CreateGaussianPosteriorMatrix(posterior,  model.parameter->t[level], center, scale); 

        if (!proposal_matrix || !posterior_matrix)
        {
                cerr << "There are no proposal or posterior draws.\n";
                abort();
        }

        int in_P1, in_P2;
        double mdd_mueller = ComputeLogMarginalDensity_Mueller(proposal_matrix, posterior_matrix, &in_P1, &in_P2);
        double mdd_bridge = ComputeLogMarginalDensity_Bridge(proposal_matrix, posterior_matrix);

	FreeMatrix(proposal_matrix); 
	FreeMatrix(posterior_matrix); 
	return mdd_bridge; 
}

double EstimateLogMDD(CEquiEnergyModel &model, int level, int proposal_type)
{
	vector<CSampleIDWeight > posterior; 
	bool unstructured = true; 
	 int data_size = model.current_sample.GetSize_Data(); 
	if (!model.storage->DrawAllSample(level, posterior, unstructured, data_size)) 
        {
                cerr << "EstimateLogMDD:: error occurred when loading all samples.\n"; 
                abort();
        }

	TDenseVector center(posterior[0].data.dim, 0.0); 
	TDenseMatrix scale(posterior[0].data.dim, posterior[0].data.dim, 0.0); 
	if(!GetCenterScaleFromSample(posterior, center, scale))
	{
		cerr << "Error occurred in GetCenterScaleFromSample()\n"; 
		abort(); 
	}
	TDenseVector R = GetRadiusFromSample(posterior, center, scale); 

	TVector center_vector = CreateVector(center.dim); 
	TMatrix scale_matrix = CreateMatrix(scale.rows, scale.cols); 
	TVector R_vector = CreateVector(R.dim); 
	for (int i=0; i<center.dim; i++)
		ElementV(center_vector, i) = center[i]; 
	for (int i=0; i<scale.rows; i++)
		for (int j=0; j<scale.cols; j++)
			ElementM(scale_matrix, i, j) = scale(i,j); 
	for (int i=0; i<R.dim; i++)
		ElementV(R_vector, i) = R[i]; 
	
	double low=0.1, high=0.9; 
	TElliptical *elliptical=(TElliptical*)NULL; 
	switch (proposal_type)
        {
                case USE_GAUSSIAN:
                        elliptical=CreateEllipticalFromPosterior_Gaussian(R_vector,center.dim,center_vector,scale_matrix);
                        break;
                case USE_POWER:
                        elliptical=CreateEllipticalFromPosterior_Power(R_vector,center.dim,center_vector,scale_matrix, high);
                        break;
                case USE_TRUNCATED_POWER:
                        elliptical=CreateEllipticalFromPosterior_TruncatedPower(R_vector,center.dim,center_vector,scale_matrix,low,high);
                        break;
                case USE_STEP:
                        elliptical=CreateEllipticalFromPosterior_Step(R_vector,center.dim,center_vector,scale_matrix,low,high);
                        break;
                case USE_TRUNCATED_GAUSSIAN:
                        elliptical=CreateEllipticalFromPosterior_TruncatedGaussian(R_vector,center.dim,center_vector,scale_matrix,low,high);
                        break;
                case USE_UNIFORM:
                        elliptical=CreateEllipticalFromPosterior_Uniform(R_vector,center.dim,center_vector,scale_matrix,low,high);
                        break;
                default:
                        cerr << "proposal type can only be 1 (Gaussian), 2(power), 3 (truncated power), 4 (step), 5 (truncated gaussian), 6 (uniform).\n";
                        abort();
                        break;
        }

	// TMatrix proposal_matrix=CreateProposalMatrix(posterior.size(), model, elliptical); 
	// TMatrix posterior_matrix=CreatePosteriorMatrix(posterior, elliptical); 
	TMatrix proposal_matrix=CreateProposalMatrix(posterior.size(), model, model.parameter->t[level], elliptical); 
	TMatrix posterior_matrix=CreatePosteriorMatrix(posterior, model.parameter->t[level], elliptical); 

	if (!proposal_matrix || !posterior_matrix)
	{
		cerr << "There are no proposal or posterior draws.\n"; 
		abort(); 
	}

	int in_P1, in_P2; 
	double mdd_mueller = ComputeLogMarginalDensity_Mueller(proposal_matrix, posterior_matrix, &in_P1, &in_P2); 
	double mdd_bridge = ComputeLogMarginalDensity_Bridge(proposal_matrix, posterior_matrix); 
	
	// cout << "logMDD mueller at level " << level << " using truncated power as proposal is " << mdd_mueller << endl; 
	// cout << "logMDD bridge at level " << level << " using truncated power as proposal is " << mdd_bridge << endl; 

	FreeVector(center_vector); 
	FreeMatrix(scale_matrix); 
	FreeVector(R_vector); 
	FreeMatrix(proposal_matrix); 
	FreeMatrix(posterior_matrix); 
	FreeElliptical(elliptical);

	return mdd_bridge; 
}


double EstimateLogMDD(CEquiEnergyModel &model, int level, int previous_level,  double logMDD_previous)
{
	vector<CSampleIDWeight> proposal, posterior; 

	bool unstructured = true; 
	int data_size = model.current_sample.GetSize_Data();
	if (!model.storage->DrawAllSample(previous_level, proposal, unstructured, data_size) || !model.storage->DrawAllSample(level, posterior, unstructured, data_size)) 
	{
		cerr << "EstimateLogMDD:: error occurred when loading all samples.\n"; 
		abort(); 
	}

	// dimension 0: proposal density
	// dimension 1: posterior kernel
	TMatrix proposal_value=CreateMatrix(proposal.size(), 2); 
	TMatrix posterior_value=CreateMatrix(posterior.size(), 2);  

	for(int i=0; i<(int)proposal.size(); i++)
	{
		ElementM(proposal_value, i, 0) = proposal[i].weight/model.parameter->t[previous_level]-logMDD_previous; 
        	ElementM(proposal_value, i, 1) = proposal[i].weight/model.parameter->t[level];
		// ElementM(proposal_value, i, 0) = proposal[i].weight-logMDD_previous; 
        	// ElementM(proposal_value, i, 1) = proposal[i].weight;
	}
	for(int i=0; i<(int)posterior.size(); i++)
	{
		ElementM(posterior_value, i, 0) = posterior[i].weight/model.parameter->t[previous_level]-logMDD_previous; 
                ElementM(posterior_value, i, 1) = posterior[i].weight/model.parameter->t[level];
		// ElementM(posterior_value, i, 0) = posterior[i].weight-logMDD_previous; 
                // ElementM(posterior_value, i, 1) = posterior[i].weight;
	}
	
	int in_P1, in_P2; 
	double logmdd = ComputeLogMarginalDensity_Bridge(proposal_value, posterior_value); 
	double logMDD = ComputeLogMarginalDensity_Mueller(proposal_value, posterior_value, &in_P1, &in_P2); 
	
	// cout << "logMDD mueller at level " << level <<  " using draws of level " << previous_level << " as proposal is " << logMDD << endl; 
	// cout << "logMDD bridge at level " << level << " using draws of previuos level " << previous_level << " as proposal is " << logmdd << endl; 
	return logMDD; 
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


