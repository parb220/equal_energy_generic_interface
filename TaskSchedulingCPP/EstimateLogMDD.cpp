#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
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

bool compare_CSampleIDWeight_BasedOnEnergy(const CSampleIDWeight &i, const CSampleIDWeight &j); 
bool compare_CSampleIDWeight_BasedOnID(const CSampleIDWeight &i, const CSampleIDWeight &j); 
double EstimateLogMDD(CEquiEnergyModel &model, int level, int previous_level, double logMDD_previous); 
double EstimateLogMDD(CEquiEnergyModel &model, int level, int proposal_type);
double EstimateLogMDD_gaussian(CEquiEnergyModel &model, int level); 

vector<double>EstimateLogMDD(CEquiEnergyModel &model, int level, int proposal_type, int nGroup); 
vector<double>EstimateLogMDD_gaussian(CEquiEnergyModel &model, int level, int nGroup);
 
double LogMDD(const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, double t, int proposal_type); 
double LogMDD_gaussian(const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, double t); 
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
	return LogMDD_gaussian(posterior, model, model.parameter->t[level]); 
}

double LogMDD_gaussian(const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, double t)
{	
	TDenseVector center(posterior[0].data.dim, 0.0);
        TDenseMatrix scale(posterior[0].data.dim, posterior[0].data.dim, 0.0);
        if(!GetCenterScaleFromSample(posterior, center, scale))
        {
                cerr << "Error occurred in GetCenterScaleFromSample()\n";
                abort();
        }

	TMatrix proposal_matrix = CreateGaussianProposalMatrix(posterior.size(), model, t, center, scale); 
	TMatrix posterior_matrix = CreateGaussianPosteriorMatrix(posterior,  t, center, scale); 

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
	return LogMDD(posterior, model, model.parameter->t[level], proposal_type); 
}

double LogMDD(const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, double t,  int proposal_type)
{
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
	TMatrix proposal_matrix=CreateProposalMatrix((int)posterior.size(), model, t, elliptical); 
	TMatrix posterior_matrix=CreatePosteriorMatrix(posterior, t, elliptical); 

	if (!proposal_matrix || !posterior_matrix)
	{
		cerr << "There are no proposal or posterior draws.\n"; 
		abort(); 
	}

	int in_P1, in_P2; 
	double mdd_mueller = ComputeLogMarginalDensity_Mueller(proposal_matrix, posterior_matrix, &in_P1, &in_P2); 
	double mdd_bridge = ComputeLogMarginalDensity_Bridge(proposal_matrix, posterior_matrix); 
	
	FreeVector(center_vector); 
	FreeMatrix(scale_matrix); 
	FreeVector(R_vector); 
	FreeMatrix(proposal_matrix); 
	FreeMatrix(posterior_matrix); 
	FreeElliptical(elliptical);

	return mdd_bridge; 
}

double LowerBoundEffectiveSampleSize(CEquiEnergyModel &model, int level, int previous_level)
{
	vector<CSampleIDWeight> proposal; 

	bool unstructured = true; 
	int data_size = model.current_sample.GetSize_Data();
	if (!model.storage->DrawAllSample(previous_level, proposal, unstructured, data_size) ) 
	{
		cerr << "LowerBoundEffectiveSampleSize():: error occurred when checking convergency.\n"; 
		abort(); 
	}

	vector<double> weight(proposal.size(), 0.0); 
	double sum_weight=0.0; 
	for(int i=0; i<(int)proposal.size(); i++)
	{
		if (previous_level == model.parameter->number_energy_level)
			weight[i] = proposal[i].weight/model.parameter->t[level]-model.StudentT_LogPDF(proposal[i]);
		else 
			weight[i] = proposal[i].weight/model.parameter->t[level]-proposal[i].weight/model.parameter->t[previous_level]; 
		if (i==0)
			sum_weight = weight[i]; 
		else
			sum_weight = AddLogs(sum_weight, weight[i]); 
	}
	sort(weight.begin(), weight.end()); 

	double LB_ESS = 0.0; 
	for (int i=0; i<(int)proposal.size(); i++)
		LB_ESS += exp(2.0*(weight[i]-sum_weight)); 
	
	return 1.0/LB_ESS; 
}

vector<double> EffectiveSampleSize(CEquiEnergyModel &model, int level)
{
	vector<CSampleIDWeight> sample; 
	bool unstructured = true;
	int data_size = model.current_sample.GetSize_Data();
	if (!model.storage->DrawAllSample(level, sample, unstructured, data_size) )
        {
                cerr << "EffectiveSampleSize():: error occurred when checking convergency.\n";
                abort();
        }
	sort(sample.begin(), sample.end(), compare_CSampleIDWeight_BasedOnID);

	vector<TDenseVector > group_mean, group_variance; 
	TDenseVector mean;
	int counter =0; 
	for (int i=0; i<(int)sample.size(); i++)
	{
		if (i == 0)
		{
			group_mean.push_back(sample[i].data); 
			group_variance.push_back(DotMultiply(sample[i].data, sample[i].data)); 
			counter = 1; 
			mean = sample[i].data; 
		}
		else if (sample[i].id > sample[i-1].id)
		{
			group_mean[group_mean.size()-1] = group_mean[group_mean.size()-1] * (1.0/(double)counter); 
			group_variance[group_variance.size()-1] = group_variance[group_variance.size()-1] * (1.0/(double)counter) - DotMultiply(group_mean[group_mean.size()-1], group_mean[group_mean.size()-1]); 
			group_mean.push_back(sample[i].data);
                        group_variance.push_back(DotMultiply(sample[i].data, sample[i].data));
			counter = 1; 
			mean = mean + sample[i].data;
		}
		else
		{
			group_mean[group_mean.size()-1] = group_mean[group_mean.size()-1] + sample[i].data; 
			group_variance[group_variance.size()-1] = group_variance[group_variance.size()-1]  + DotMultiply(sample[i].data, sample[i].data);
			counter ++; 
			mean = mean + sample[i].data;
		}
	}
	group_mean[group_mean.size()-1] = group_mean[group_mean.size()-1] * (1.0/(double)counter);
        group_variance[group_variance.size()-1] = group_variance[group_variance.size()-1] * (1.0/(double)counter) - DotMultiply(group_mean[group_mean.size()-1], group_mean[group_mean.size()-1]);

	mean = mean * (1.0/(double)sample.size()); 

	TDenseVector B(mean.dim,0.0), W(mean.dim,0.0); 
	for (int i=0; i<(int)group_mean.size(); i++)
	{
		B = B + DotMultiply(group_mean[i]-mean, group_mean[i]-mean); 
		W = W + group_variance[i]; 
	}	
	B = B * ((double)counter/(double)group_mean.size()); 
	W = W * (1.0/(double)group_mean.size()); 

	vector<double> ESS(B.dim, 0.0); 
	for (int i=0; i<B.dim; i++)
		ESS[i] = (double)sample.size()/(B[i]/W[i]); 

	return ESS; 
}

double CheckConvergency (CEquiEnergyModel &model, int level, int previous_level,  double convergency_previous, double &average_consistency, double &std_consistency, double &LB_ESS)
{
	vector<CSampleIDWeight> proposal; 

	bool unstructured = true; 
	int data_size = model.current_sample.GetSize_Data();
	if (!model.storage->DrawAllSample(previous_level, proposal, unstructured, data_size) ) 
	{
		cerr << "CheckConvergency:: error occurred when checking convergency.\n"; 
		abort(); 
	}

	sort(proposal.begin(), proposal.end(), compare_CSampleIDWeight_BasedOnID);
	vector<double> weight(proposal.size(), 0.0); 
	for(int i=0; i<(int)proposal.size(); i++)
	{
		if (previous_level == model.parameter->number_energy_level)
			weight[i] = proposal[i].weight/model.parameter->t[level]-model.StudentT_LogPDF(proposal[i]);
		else 
			weight[i] = proposal[i].weight/model.parameter->t[level]-proposal[i].weight/model.parameter->t[previous_level]; 
	}

	vector<double>group_consistency; 
	double consistency, sum_weight; 
	int counter =0; 
	for (int i=0; i<(int)(proposal.size()); i++)
	{
		if (i==0)
		{ 
			group_consistency.push_back(weight[i]); 
			counter = 1; 
			consistency = weight[i]; 
		}
		else if (proposal[i].id > proposal[i-1].id)
		{
			group_consistency[group_consistency.size()-1] = group_consistency[group_consistency.size()-1] + convergency_previous - log((double)counter); 
			group_consistency.push_back(weight[i]); 
			counter = 1; 
			consistency = AddLogs(consistency, weight[i]); 
		}
		else 
		{
			group_consistency[group_consistency.size()-1] = AddLogs(group_consistency[group_consistency.size()-1], weight[i]);  
			counter ++; 
			consistency = AddLogs(consistency, weight[i]);
		}
		if ( i==0)
			sum_weight = weight[i]; 
		else 
			sum_weight = AddLogs(sum_weight, weight[i]);
	}
	group_consistency[group_consistency.size()-1] = group_consistency[group_consistency.size()-1] + convergency_previous- log((double)counter);
	consistency = consistency + convergency_previous- log((double)(proposal.size())); 
	
	average_consistency=0.0; 
	std_consistency=0.0; 
	for (int i=0; i<(int)(group_consistency.size()); i++)
	{
		average_consistency += group_consistency[i]; 
		std_consistency += group_consistency[i] * group_consistency[i]; 
	}	
	average_consistency = average_consistency/(double)(group_consistency.size()); 
	std_consistency = sqrt(std_consistency/(double)(group_consistency.size())-average_consistency*average_consistency); 
	
	LB_ESS = 0.0; 
	for (int i=0; i<(int)proposal.size(); i++)
		LB_ESS += exp(2.0*(weight[i]-sum_weight)); 
	LB_ESS = 1.0/LB_ESS;  

	return consistency; 
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

double CheckLogMDDConvergency(CEquiEnergyModel &model, int level, double t, int proposal_type, double &average_logMDD, double &std_logMDD)
{
	vector<CSampleIDWeight> sample;
        bool unstructured = true;
        int data_size = model.current_sample.GetSize_Data();
        if (!model.storage->DrawAllSample(level, sample, unstructured, data_size) )
        {
                cerr << "CheckConvergency:: error occurred when checking convergency.\n";
                abort();
        }
        sort(sample.begin(), sample.end(), compare_CSampleIDWeight_BasedOnID);

	vector<double>group_logMDD;
        average_logMDD = 0.0, std_logMDD = 0.0; 
        int start_index, counter =0;
        for (int i=0; i<(int)(sample.size()); i++)
        {
                if (i==0)
		{
			start_index = 0; 
                        counter = 1;
		}
                else if (sample[i].id > sample[i-1].id)
                {
			group_logMDD.push_back(LogMDD(vector<CSampleIDWeight>(sample.begin()+start_index, sample.begin()+start_index+counter), model, t, proposal_type));
			start_index = counter; 
			counter = 1; 
                }
                else
                        counter ++;
        }
	group_logMDD.push_back(LogMDD(vector<CSampleIDWeight>(sample.begin()+start_index, sample.begin()+start_index+counter), model, t, proposal_type)); 

	/*sort(group_logMDD.begin(), group_logMDD.end()); 
	cout << "logMDD per j" ; 
	for (int i=0; i<(int)(group_logMDD.size()); i++)
		cout << "\t" << setprecision(20) << group_logMDD[i]; 
	cout << endl; */
        for (int i=0; i<(int)(group_logMDD.size()); i++)
        {
                average_logMDD += group_logMDD[i];
                std_logMDD += group_logMDD[i] * group_logMDD[i];
        }
        average_logMDD = average_logMDD/(double)(group_logMDD.size());
        std_logMDD = sqrt(std_logMDD/(double)(group_logMDD.size())-average_logMDD*average_logMDD);
        return LogMDD(sample, model, t, proposal_type); 
}

vector<double>EstimateLogMDD(CEquiEnergyModel &model, int level, int proposal_type, int nGroup)
{
	vector<CSampleIDWeight> posterior; 
	bool unstructured = true;
        int data_size = model.current_sample.GetSize_Data();
	if(!model.storage->DrawAllSample(level, posterior, unstructured, data_size))
        {
                cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
                abort();
        }
	sort(posterior.begin(), posterior.end(), compare_CSampleIDWeight_BasedOnEnergy);
	random_shuffle(posterior.begin(), posterior.end()); 

	int group_size = ceil((double)posterior.size()/(double)nGroup); 
	int begin_index = 0, end_index = group_size; 
	vector<CSampleIDWeight> group_sample; 
	vector<double> logMDD(nGroup, 0.0); 
	for (int i=0; i<nGroup; i++)
	{
		group_sample.clear(); 
		for (int j=begin_index; j<end_index; j++)
			group_sample.push_back(posterior[j]); 
		logMDD[i] = LogMDD(group_sample, model, model.parameter->t[level], proposal_type); 
		begin_index = end_index; 
		end_index = end_index + group_size < (int)(posterior.size()) ? end_index + group_size : (int)(posterior.size()); 
	}

	return logMDD; 
}

vector<double>EstimateLogMDD_gaussian(CEquiEnergyModel &model, int level, int nGroup)
{
	vector<CSampleIDWeight> posterior;
        bool unstructured = true;
        int data_size = model.current_sample.GetSize_Data();
        if(!model.storage->DrawAllSample(level, posterior, unstructured, data_size))
        {
                cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
                abort();
        }
        sort(posterior.begin(), posterior.end(), compare_CSampleIDWeight_BasedOnEnergy);
	random_shuffle(posterior.begin(), posterior.end()); 

        int group_size = ceil((double)posterior.size()/(double)nGroup);
        int begin_index = 0, end_index = group_size;
        vector<CSampleIDWeight> group_sample;
        vector<double> logMDD(nGroup, 0.0);
        for (int i=0; i<nGroup; i++)
        {
                group_sample.clear();
                for (int j=begin_index; j<end_index; j++)
                        group_sample.push_back(posterior[j]);
                logMDD[i] = LogMDD_gaussian(group_sample, model, model.parameter->t[level]);
                begin_index = end_index;
                end_index = end_index + group_size < (int)(posterior.size()) ? end_index + group_size : (int)(posterior.size());
        }
        return logMDD;
}

