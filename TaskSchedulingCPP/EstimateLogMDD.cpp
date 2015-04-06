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
#include "EstimateLogMDD.hpp"

using namespace std; 

bool compare_CSampleIDWeight_BasedOnEnergy(const CSampleIDWeight &i, const CSampleIDWeight &j); 
bool compare_CSampleIDWeight_BasedOnID(const CSampleIDWeight &i, const CSampleIDWeight &j); 

// Function used in LogMDD calculationg (no need for external interface)
TMatrix CreateProposalMatrix(int ndraws, CEquiEnergyModel &model, double t, const TElliptical *elliptical, int posterior_type=POSTERIOR_HEATED); 
TMatrix CreatePosteriorMatrix(const vector<CSampleIDWeight> &sample, double t, const TElliptical *elliptical, int posterior_type=POSTERIOR_HEATED);
TMatrix CreateGaussianProposalMatrix(int ndraws, CEquiEnergyModel &model, double t, const TDenseVector &center, const TDenseMatrix &scale, int posterior_type=POSTERIOR_HEATED);
TMatrix CreateGaussianPosteriorMatrix(const vector<CSampleIDWeight> &sample, double t, const TDenseVector &center, const TDenseMatrix &scale, int posterior_type=POSTERIOR_HEATED);

TDenseVector GetRadiusFromSample(const vector<CSampleIDWeight> &sample, const TDenseVector &center, const TDenseMatrix &scale); 
bool GetCenterScaleFromSample(const vector<CSampleIDWeight> &sample, TDenseVector &center, TDenseMatrix &scale, double t, int posterior_type=POSTERIOR_HEATED); 


TMatrix CreateGaussianPosteriorMatrix(const vector<CSampleIDWeight> &sample, double t, const TDenseVector &center, const TDenseMatrix &scale, int posterior_type)
{
	TMatrix posterior = CreateMatrix(sample.size(), 2); 
	double logdensity; 
	TDenseVector y(sample[0].data.dim, 0.0); 
	TDenseMatrix scale_inverse = Inverse(scale, SOLVE_SVD); 
	
	for (int i=0; i<(int)sample.size(); i++)
	{
		y = scale_inverse*(sample[i].data-center);  // y = Inverse(scale)*(x-center); 
		logdensity = -0.918938533204673 -0.5*y[0]*y[0];
		for (int j=1; j<y.dim; j++)
                        logdensity += -0.918938533204673 -0.5*y[j]*y[j];
		ElementM(posterior,i,0)=logdensity;
		if (posterior_type == POSTERIOR_HEATED)
			ElementM(posterior, i, 1) = sample[i].weight/t; 
		else if (posterior_type == LIKELIHOOD_HEATED)
			ElementM(posterior, i, 1) = sample[i].reserved/t + sample[i].weight-sample[i].reserved; 
		else if (posterior_type == PRIOR_ONLY)
			ElementM(posterior, i, 1) = sample[i].weight-sample[i].reserved; 
	}
	return posterior; 
}

TMatrix CreateGaussianProposalMatrix(int ndraws, CEquiEnergyModel &model, double t, const TDenseVector &center, const TDenseMatrix &scale, int posterior_type)
{
	if (ndraws <= 0)
		return (TMatrix)NULL; 
	
	TMatrix proposal = CreateMatrix(ndraws, 2); 
	CSampleIDWeight sample; 

	TDenseVector x(center.dim,0.0), y(center.dim, 0.0); 
	double log_density; 
	for (int i=0; i<ndraws; i++)
	{
		sample.data = RandomNormalVector(center.dim); 
		log_density = 0; 
		for (int j=0; j<center.dim; j++)
			log_density += -0.918938533204673 -0.5*sample.data[j]*sample.data[j];

		sample.data = scale * sample.data + center; // y = center + scale *x; 
		sample.DataChanged();
		model.log_posterior_function(sample);  
		ElementM(proposal, i, 0) = log_density;
		if (posterior_type == POSTERIOR_HEATED)
                	ElementM(proposal, i, 1) = sample.weight/t; 
		else if (posterior_type == LIKELIHOOD_HEATED)
			ElementM(proposal, i, 1) = sample.reserved/t+(sample.weight-sample.reserved); 
		else if (posterior_type == PRIOR_ONLY)
			ElementM(proposal, i, 1) = sample.weight-sample.reserved; 
	}
	return proposal; 
}


TMatrix CreateProposalMatrix(int ndraws, CEquiEnergyModel &model, double t, const TElliptical *elliptical, int posterior_type)
{
	if (ndraws <= 0)
		return (TMatrix)NULL; 

	TMatrix proposal = CreateMatrix(ndraws, 2); 
	int dim = elliptical ? elliptical->dim : 0;
	double *draw = dim>0 ? new double[dim] : NULL; 

	CSampleIDWeight sample; 
	sample.data.Resize(dim); 
	double log_density, log_posterior;
	for (int i=0; i<ndraws; i++)
	{
		log_density = draw ?  LogDensityElliptical_Radius(DrawElliptical(draw, (TElliptical *)elliptical), (TElliptical *)elliptical) : 0.0;
		ElementM(proposal, i, 0) = log_density; 
	
		for (int j=0; j<dim; j++)
			sample.data[j] = draw[j]; 
		sample.DataChanged(); 
		model.log_posterior_function(sample); 

		if (posterior_type == POSTERIOR_HEATED)
			log_posterior = sample.weight/t; 
		else if (posterior_type == LIKELIHOOD_HEATED)
			log_posterior = sample.reserved/t + (sample.weight-sample.reserved); 
		else if (posterior_type == PRIOR_ONLY)
			log_posterior = sample.weight-sample.reserved; 
		ElementM(proposal, i, 1) = log_posterior; 
	}

	if (draw)
		delete [] draw; 
	return proposal; 
}

TMatrix CreatePosteriorMatrix(const vector<CSampleIDWeight> &sample, double t, const TElliptical *elliptical, int posterior_type)
{
	TMatrix posterior=CreateMatrix(sample.size(),2); 
	int dim = elliptical ? elliptical->dim:0; 

	double log_density; 
	for (int i=0; i<(int)sample.size(); i++)
	{
		log_density = dim> 0 ? LogDensityElliptical_Draw(sample[i].data.vector, (TElliptical *)elliptical) : 0.0; 
		ElementM(posterior, i, 0) = log_density; 
		if (posterior_type == POSTERIOR_HEATED)
			ElementM(posterior, i, 1) = sample[i].weight/t;
		else if (posterior_type == LIKELIHOOD_HEATED)
			ElementM(posterior, i, 1) = sample[i].reserved/t + (sample[i].weight-sample[i].reserved); 
		else if (posterior_type == PRIOR_ONLY)
			ElementM(posterior, i, 1) = sample[i].weight - sample[i].reserved; 
	}
	return posterior; 
}


TDenseVector GetRadiusFromSample(const vector<CSampleIDWeight> &sample, const TDenseVector &center, const TDenseMatrix &scale)
{
	TDenseVector Radii(sample.size(), 0.0); 
	TDenseMatrix quadratic_form = MultiplyTranspose(scale, scale); 
	quadratic_form.Inverse(SOLVE_SVD);
	quadratic_form = (quadratic_form + Transpose(quadratic_form)) * 0.5; // force symmetry
	
	TDenseVector theta; 
	for (int i=0; i<(int)sample.size(); i++)
	{
		theta = sample[i].data - center; 
		Radii[i] =  sqrt(InnerProduct(theta, theta, quadratic_form)); 
	}
	return Radii; 
}

bool GetCenterScaleFromSample(const vector<CSampleIDWeight> &sample, TDenseVector &center, TDenseMatrix &scale, double t, int posterior_type)
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
		if (i==0 || (posterior_type == POSTERIOR_HEATED && sample[i].weight/t > max) || (posterior_type == LIKELIHOOD_HEATED && sample[i].reserved/t+(sample[i].weight-sample[i].reserved) > max) || (posterior_type == PRIOR_ONLY && (sample[i].weight-sample[i].reserved) > max))
		{
			if (posterior_type == POSTERIOR_HEATED)
				max = sample[i].weight/t; 
			else if (posterior_type == LIKELIHOOD_HEATED)
				max = sample[i].reserved/t + (sample[i].weight-sample[i].reserved); 
			else if (posterior_type == PRIOR_ONLY) 
				max = sample[i].weight-sample[i].reserved; 
			mode = sample[i];  
		}
		sample_sum = sample_sum + sample[i].data; 
		sample_square = sample_square + OuterProduct(sample[i].data, sample[i].data); 
	}
	
	sample_sum = sample_sum * (1.0/(double)sample.size()); 
	sample_square = sample_square * (1.0/(double)sample.size()); 	

	// using mode
	center = mode.data; 
	sample_square = sample_square + OuterProduct(mode.data, mode.data); 
	sample_square = sample_square - OuterProduct(sample_sum, mode.data) - OuterProduct(mode.data, sample_sum);
	// using center
	// sample_square = sample_square - OuterProduct(sample_sum, sample_sum); 	
 
	sample_square = (sample_square+Transpose(sample_square))*0.5; 
	
	/*scale = Cholesky(sample_square); 
	scale = Transpose(scale); */
	
	// Eig analysis
	TDenseVector eValue(sample_square.rows,0.0); 
	TDenseMatrix eVector(sample_square.rows, sample_square.cols, 0.0); 
	Eig(eValue, eVector, sample_square); 
	for (int i=0; i<eValue.dim; i++)
		eValue[i] = sqrt(eValue[i]); 
	scale = MultiplyTranspose(eVector * DiagonalMatrix(eValue), eVector); 
	return true; ; 
}

double EstimateLogMDD_gaussian(CEquiEnergyModel &model, int stage, int posterior_type)
{
	vector<CSampleIDWeight > posterior;
        bool unstructured = true;
        int data_size = model.current_sample.GetSize_Data();
        if (!model.storage->DrawAllSample(stage, posterior, unstructured, data_size))
        {
                cerr << "EstimateLogMDD:: error occurred when loading all samples.\n";
                abort();
        }
	double t = model.parameter->t[stage]; 
	return LogMDD_gaussian(posterior, model, t, posterior_type);
}

double LogMDD_gaussian(const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, double t, int posterior_type)
{ 
	TDenseVector center(posterior[0].data.dim, 0.0);
        TDenseMatrix scale(posterior[0].data.dim, posterior[0].data.dim, 0.0);
        if(!GetCenterScaleFromSample(posterior, center, scale, t, posterior_type))
        {
                cerr << "Error occurred in GetCenterScaleFromSample()\n";
                abort();
        }

	TMatrix proposal_matrix = CreateGaussianProposalMatrix(posterior.size(), model, t, center, scale, posterior_type); 
	TMatrix posterior_matrix = CreateGaussianPosteriorMatrix(posterior,  t, center, scale, posterior_type); 

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


vector<double> WeightFromGaussianSample(int N, const TDenseVector &center, const TDenseMatrix &scale, CEquiEnergyModel &model, double temperature, int posterior_type)
{
	vector<double> weight(N,0.0); 
	CSampleIDWeight sample; 
	sample.data.Resize(center.dim); 
	for (int i=0; i<N; i++)
	{
		sample.data = RandomNormalVector(center.dim);
		sample.data = scale * sample.data + center; 
		sample.DataChanged(); 
		model.log_prior_function(sample); 
	
		if (posterior_type == POSTERIOR_HEATED)
			weight[i] = sample.weight/temperature; 
		else if (posterior_type == LIKELIHOOD_HEATED)
			weight[i] = sample.reserved/temperature + sample.weight-sample.reserved;  
		else if (posterior_type == PRIOR_ONLY)
			weight[i] = sample.weight-sample.reserved; 
		weight[i] += 0.918938533204673*(double)center.dim; // 0.5*log(2*pi) = 0.918938533204673
		for (int j=0; j<center.dim; j++)
			weight[i] += 0.5*sample.data[j]*sample.data[j];  
	}
	return weight; 
}

double EstimateLogMDD(CEquiEnergyModel &model, int stage, int proposal_type, int posterior_type)
{
	vector<CSampleIDWeight > posterior; 
	bool unstructured = true; 
	 int data_size = model.current_sample.GetSize_Data(); 
	if (!model.storage->DrawAllSample(stage, posterior, unstructured, data_size)) 
        {
                cerr << "EstimateLogMDD:: error occurred when loading all samples.\n"; 
                abort();
        }
	double t = model.parameter->t[stage]; 
	return LogMDD(posterior, model, t, proposal_type, posterior_type);
}

double LogMDD(const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, double t, int proposal_type, int posterior_type)
{	
	TDenseVector center(posterior[0].data.dim, 0.0);
        TDenseMatrix scale(posterior[0].data.dim, posterior[0].data.dim, 0.0);
        if(!GetCenterScaleFromSample(posterior, center, scale, t, posterior_type))
        {
                cerr << "Error occurred in GetCenterScaleFromSample()\n";
                abort();
        } 
	 
	if (proposal_type == WEIGHTED_WITH_GAUSIAN_SAMPLES)
	{
		vector<double> weight = WeightFromGaussianSample((int)posterior.size(), center, scale, model, t, posterior_type);
		double sum_weight = weight[0]; 
		for (int i=1; i<(int)weight.size(); i++)
                        sum_weight = AddLogs(sum_weight, weight[i]);
		return sum_weight; 
	}	
	else 
	{
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

		TMatrix proposal_matrix=CreateProposalMatrix((int)posterior.size(), model, t, elliptical, posterior_type); 
		TMatrix posterior_matrix=CreatePosteriorMatrix(posterior, t, elliptical, posterior_type); 

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
}

double LowerBoundEffectiveSampleSize(CEquiEnergyModel &model, int stage, int previous_stage, int posterior_type)
{
	vector<CSampleIDWeight> proposal; 

	bool unstructured = true; 
	int data_size = model.current_sample.GetSize_Data();
	if (!model.storage->DrawAllSample(previous_stage, proposal, unstructured, data_size) ) 
	{
		cerr << "LowerBoundEffectiveSampleSize():: error occurred when checking convergency.\n"; 
		abort(); 
	}

	vector<double> weight(proposal.size(), 0.0); 
	double sum_weight=0.0; 
	for(int i=0; i<(int)proposal.size(); i++)
	{
		if (previous_stage == model.parameter->number_energy_stage)
		{
			if (posterior_type == POSTERIOR_HEATED)
				weight[i] = proposal[i].weight/model.parameter->t[stage]-(proposal[i].weight-proposal[i].reserved); 
			else if (posterior_type = LIKELIHOOD_HEATED)
				weight[i] = proposal[i].reserved/model.parameter->t[stage]; 
			else if (posterior_type == PRIOR_ONLY)
				weight[i]  = 0; 
		}
		else 
		{
			if (posterior_type == POSTERIOR_HEATED)
				weight[i] = proposal[i].weight/model.parameter->t[stage]-proposal[i].weight/model.parameter->t[previous_stage]; 
			else if (posterior_type = LIKELIHOOD_HEATED)
				weight[i] = proposal[i].reserved/model.parameter->t[stage]-proposal[i].reserved/model.parameter->t[previous_stage]; 
			else if (posterior_type == PRIOR_ONLY)
				weight[i]  = 0; 
		}			

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

vector<double> EffectiveSampleSize(vector<CSampleIDWeight> &sample) // CEquiEnergyModel &model, int stage, bool if_normalization)
{
	/*vector<CSampleIDWeight> sample; 
	bool unstructured = true;
	int data_size = model.current_sample.GetSize_Data();
	if (!model.storage->DrawAllSample(stage, sample, unstructured, data_size) )
        {
                cerr << "EffectiveSampleSize():: error occurred when checking convergency.\n";
                abort();
        }*/
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

double CheckConvergency (CEquiEnergyModel &model, int stage, int previous_stage,  double convergency_previous, double &average_consistency, double &std_consistency, double &LB_ESS, int posterior_type)
{
	vector<CSampleIDWeight> proposal; 

	bool unstructured = true; 
	int data_size = model.current_sample.GetSize_Data();
	if (!model.storage->DrawAllSample(previous_stage, proposal, unstructured, data_size) ) 
	{
		cerr << "CheckConvergency:: error occurred when checking convergency.\n"; 
		abort(); 
	}

	sort(proposal.begin(), proposal.end(), compare_CSampleIDWeight_BasedOnID);
	vector<double> weight(proposal.size(), 0.0); 
	for(int i=0; i<(int)proposal.size(); i++)
	{
		if (previous_stage == model.parameter->number_energy_stage)
		{
			if (posterior_type == POSTERIOR_HEATED)
				weight[i] = proposal[i].weight/model.parameter->t[stage] - (proposal[i].weight - proposal[i].reserved); 
			else if (posterior_type == LIKELIHOOD_HEATED)
				weight[i] = proposal[i].reserved/model.parameter->t[stage]; 
			else if (posterior_type == PRIOR_ONLY)
				weight[i] = 0.0; 
		}
		else 
		{
			if (posterior_type == POSTERIOR_HEATED)
				weight[i] = proposal[i].weight/model.parameter->t[stage] - proposal[i].weight/model.parameter->t[previous_stage]; 
			else if (posterior_type == LIKELIHOOD_HEATED)
				weight[i] = proposal[i].reserved/model.parameter->t[stage] - proposal[i].reserved/model.parameter->t[previous_stage]; 
			else if (posterior_type == PRIOR_ONLY)
				weight[i] = 0.0; 
		}
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


double EstimateLogMDD(CEquiEnergyModel &model, int stage, int previous_stage,  double logMDD_previous, int posterior_type)
{
	vector<CSampleIDWeight> proposal, posterior; 

	bool unstructured = true; 
	int data_size = model.current_sample.GetSize_Data();
	if (!model.storage->DrawAllSample(previous_stage, proposal, unstructured, data_size) || !model.storage->DrawAllSample(stage, posterior, unstructured, data_size)) 
	{
		cerr << "EstimateLogMDD:: error occurred when loading all samples.\n"; 
		abort(); 
	}
	return LogMDD(proposal, posterior, model, previous_stage, stage, logMDD_previous, posterior_type); 
}

double LogMDD(const vector<CSampleIDWeight> &proposal, const vector<CSampleIDWeight> &posterior, CEquiEnergyModel &model, int previous_stage, int current_stage, double logMDD_previous, int posterior_type)
{
	double t_previous = model.parameter->t[previous_stage], t_current = model.parameter->t[current_stage]; 

	// dimension 0: proposal density
	// dimension 1: posterior kernel
	TMatrix proposal_value=CreateMatrix(proposal.size(), 2); 
	TMatrix posterior_value=CreateMatrix(posterior.size(), 2);  

	for(int i=0; i<(int)proposal.size(); i++)
	{
		if (previous_stage == model.parameter->number_energy_stage || posterior_type == PRIOR_ONLY)
			ElementM(proposal_value, i, 0) = (proposal[i].weight - proposal[i].reserved) - logMDD_previous; 
		else if (posterior_type == POSTERIOR_HEATED)
			ElementM(proposal_value, i, 0) = proposal[i].weight/t_previous - logMDD_previous; 
		else if (posterior_type == LIKELIHOOD_HEATED)
			ElementM(proposal_value, i, 0) = (proposal[i].reserved/t_previous+proposal[i].weight-proposal[i].reserved)-logMDD_previous; 

		if (posterior_type == POSTERIOR_HEATED)
			ElementM(proposal_value, i, 1) = proposal[i].weight/t_current; 
		else if (posterior_type = LIKELIHOOD_HEATED)
        		ElementM(proposal_value, i, 1) = proposal[i].reserved/t_current+proposal[i].weight-proposal[i].reserved;
		else if (posterior_type == PRIOR_ONLY)
			ElementM(proposal_value, i, 1) = proposal[i].weight-proposal[i].reserved; 
	}
	for(int i=0; i<(int)posterior.size(); i++)
	{
		if (previous_stage == model.parameter->number_energy_stage || posterior_type == PRIOR_ONLY)
			ElementM(posterior_value, i, 0) = (posterior[i].weight-posterior[i].reserved)-logMDD_previous;
		else if (posterior_type == POSTERIOR_HEATED)
			ElementM(posterior_value, i, 0) = posterior[i].weight/t_previous - logMDD_previous; 
		else if (posterior_type == LIKELIHOOD_HEATED) 
			ElementM(posterior_value, i, 0) = (posterior[i].reserved/t_previous+posterior[i].weight-posterior[i].reserved)-logMDD_previous; 

		if (posterior_type == POSTERIOR_HEATED)
			ElementM(posterior_value, i, 1) = posterior[i].weight/t_current; 
		else if (posterior_type == LIKELIHOOD_HEATED)
                	ElementM(posterior_value, i, 1) = posterior[i].reserved/t_current+posterior[i].weight-posterior[i].reserved; 
		else if (posterior_type == PRIOR_ONLY)
			ElementM(posterior_value, i, 1) = posterior[i].weight-posterior[i].reserved; 
	}
	
	int in_P1, in_P2; 
	double logmdd = ComputeLogMarginalDensity_Bridge(proposal_value, posterior_value); 
	double logMDD = ComputeLogMarginalDensity_Mueller(proposal_value, posterior_value, &in_P1, &in_P2); 
	
	// cout << "logMDD mueller at stage " << stage <<  " using draws of stage " << previous_stage << " as proposal is " << logMDD << endl; 
	// cout << "logMDD bridge at stage " << stage << " using draws of previuos stage " << previous_stage << " as proposal is " << logmdd << endl; 
	return logMDD; 
}

double LogMDD_Importance(const vector<CSampleIDWeight> &proposal, CEquiEnergyModel &model, int previous_stage, int current_stage, double logMDD_previous, int posterior_type)
{
	double t_previous = model.parameter->t[previous_stage]; 
	double t_current = model.parameter->t[current_stage]; 

	vector<double> weight(proposal.size(),0); 
	for (int i=0; i<(int)proposal.size(); i++)
	{
		if (posterior_type == POSTERIOR_HEATED)
			weight[i] = proposal[i].weight/t_current; 
		else if (posterior_type == LIKELIHOOD_HEATED)
			weight[i] = proposal[i].reserved/t_current + (proposal[i].weight-proposal[i].reserved); 
		else if (posterior_type == PRIOR_ONLY)
			weight[i]  = proposal[i].weight-proposal[i].reserved; 

		if (previous_stage == model.parameter->number_energy_stage || posterior_type == PRIOR_ONLY)
			weight[i] -= (proposal[i].weight-proposal[i].reserved) - logMDD_previous; 
		else if (posterior_type == POSTERIOR_HEATED)
			weight[i] -= proposal[i].weight/t_previous - logMDD_previous; 
		else if (posterior_type == LIKELIHOOD_HEATED)
			weight[i] -= proposal[i].reserved/t_previous + (proposal[i].weight-proposal[i].reserved) - logMDD_previous;  
	}
	double sum_weight = weight[0];
	for (int i=1; i<(int)weight.size(); i++)
		sum_weight = AddLogs(sum_weight, weight[i]); 
	// changed by DW - must "divide" by number of observations
	return sum_weight- log(weight.size()); 
} 

double CheckLogMDDConvergency(vector<CSampleIDWeight> &sample, CEquiEnergyModel &model, int stage, double t, int proposal_type, double &average_logMDD, double &std_logMDD, int posterior_type)
{
	/*vector<CSampleIDWeight> sample;
        bool unstructured = true;
        int data_size = model.current_sample.GetSize_Data();
        if (!model.storage->DrawAllSample(stage, sample, unstructured, data_size) )
        {
                cerr << "CheckConvergency:: error occurred when checking convergency.\n";
                abort();
        }*/
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
			group_logMDD.push_back(LogMDD(vector<CSampleIDWeight>(sample.begin()+start_index, sample.begin()+start_index+counter), model, t, proposal_type, posterior_type));
			start_index = counter; 
			counter = 1; 
                }
                else
                        counter ++;
        }
	group_logMDD.push_back(LogMDD(vector<CSampleIDWeight>(sample.begin()+start_index, sample.begin()+start_index+counter), model, t, proposal_type, posterior_type)); 

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
        return LogMDD(sample, model, t, proposal_type, posterior_type); 
}

vector<double>EstimateLogMDD(CEquiEnergyModel &model, int stage, int proposal_type, int nGroup, int posterior_type)
{
	vector<CSampleIDWeight> posterior; 
	bool unstructured = true;
        int data_size = model.current_sample.GetSize_Data();
	if(!model.storage->DrawAllSample(stage, posterior, unstructured, data_size))
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
		logMDD[i] = LogMDD(group_sample, model, model.parameter->t[stage], proposal_type, posterior_type); 
		begin_index = end_index; 
		end_index = end_index + group_size < (int)(posterior.size()) ? end_index + group_size : (int)(posterior.size()); 
	}

	return logMDD; 
}

vector<double>EstimateLogMDD_gaussian(CEquiEnergyModel &model, int stage, int nGroup, int posterior_type)
{
	vector<CSampleIDWeight> posterior;
        bool unstructured = true;
        int data_size = model.current_sample.GetSize_Data();
        if(!model.storage->DrawAllSample(stage, posterior, unstructured, data_size))
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
                logMDD[i] = LogMDD_gaussian(group_sample, model, model.parameter->t[stage], posterior_type);
                begin_index = end_index;
                end_index = end_index + group_size < (int)(posterior.size()) ? end_index + group_size : (int)(posterior.size());
        }
        return logMDD;
}

