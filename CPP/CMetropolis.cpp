#include "CMetropolis.h"
#include <cmath>
#include "dw_rand.h"

using namespace std; 

void CMetropolis::BlockAdaptive(const TDenseVector &adaptive_start_point, const vector<TDenseMatrix> &B, double mid, int period, int max_period)
// Makes block-wise random-walk Metropolis draws, adjusting the block scale until the target acceptance rate 
// is hit.
//
// 	adaptive_start_point:	n*1 column vector, initial value
// 	B:			k matrices, each of which is n*b(i) consisting of the direction vectors of the i-th block
// 	mid:			target acceptance ratio
// 	period:			initial number of draws to make before recomputing scale. Increases by a factor of 2
// 				when target ratio has been reached
// 	max_period:		final number of draws to make before recomputing scale. Should be a power of 2 times
// 				period
//
// For each draw, the block i jumping kernel is as follows
// 	(1) draw z from the b(i)-dimensional standard normal distribution
// 	(2) let x = scale(i)*B{i}*z
// 	(3) if y is the current point, then y+x is the proposal point
//
// Usually b(1)+...+b(k) = n and the matrix[B{1} ... B{k}] has orthogonal columsn, though this is not checked	
{
	size_t k = B.size(); 	// number of blocks

	int n_draws =0;
	vector<int> n_accepted(k,0); 
	vector<int> begin_draws(k,0); 
	vector<int> begin_jumps(k,0); 
	vector<int> periods(k,period); 
	vector<int> end_draws(k,period); 
	
	TDenseVector previous_ratio(k,0.0); 
	TDenseVector scale(k,1.0); 
	TDenseVector best_scale(k,1.0); 
	TDenseVector low_scale(k,-1.0); 
	TDenseVector low_jump_ratio(k,-1.0); 
	TDenseVector high_scale(k,-1.0); 
	TDenseVector high_jump_ratio(k,-1.0); 

	double log_mid = log(mid); 
	double lower_bound = exp(log_mid/0.2); 	
	double upper_bound = exp(log_mid/5.0);

	vector<size_t> b(k,0); // b(i): size of the i-th block
	for (unsigned int i=0; i<k; i++)
		b[i] = B[i].cols;  

	// Beginning adaptive burn-in 
	TDenseVector x, y(adaptive_start_point.dim); 
	y.CopyContent(adaptive_start_point); 
	double log_previous = model->LogPosterior(y), log_current, new_scale, diff;
	bool done = false; 
	int check = period; 
	while (!done)
	{
		// draw metropolis blocks
		for (unsigned int i=0; i<k; i++)
		{
			x = y + scale[i]*(B[i]*RandomNormalVector(b[i]));
			log_current = model->LogPosterior(x); 

			if (log_current-log_previous >= log(dw_uniform_rnd()) )
			{
				y= x; 
				log_previous = log_current; 
				n_accepted[i] ++; 
			}
		}
		n_draws ++; 

		// rescale
		if (n_draws == check)
		{
			done = true; 
			for (unsigned int i=0; i<k; i++)
			{
				// jump ratio
				previous_ratio[i] = (double)(n_accepted[i]-begin_jumps[i])/(double)(n_draws-begin_draws[i]); 
				// recompute scale?
				if (end_draws[i] <= n_draws)
				{
					// set new low or high bounds
					if (previous_ratio[i] < mid)
					{
						low_scale[i] = scale[i]; 
						low_jump_ratio[i] = previous_ratio[i]; 
					}
					else 
					{
						high_scale[i] = scale[i]; 
						high_jump_ratio[i] = previous_ratio[i]; 
					}

					// new scale and best scale
					if (low_jump_ratio[i] < 0.0)
					{
						best_scale[i] = scale[i]; 
						if (previous_ratio[i] > upper_bound)
							new_scale = 5.0*high_scale[i]; 
						else 
							new_scale = (log_mid/log(previous_ratio(i)))*high_scale[i];
					}
					else if (high_jump_ratio[i] < 0.0)
					{
						best_scale[i] = scale[i]; 
						if (previous_ratio[i] < lower_bound)
							new_scale = 0.2*low_scale[i]; 
						else 
							new_scale=(log_mid/log(previous_ratio(i)))*low_scale[i];
					}
					else 
					{
						diff = high_jump_ratio[i]-low_jump_ratio[i]; 
						if (diff > 1.0e-6)
							new_scale=((mid-low_jump_ratio[i])*low_scale[i] + (high_jump_ratio[i]-mid)*high_scale[i])/diff;
						else 
							new_scale=0.5*(low_scale[i] + high_scale[i]);
						best_scale[i] = new_scale; 
						periods[i] = 2* periods[i]; 
						low_jump_ratio[i] = -1.0; 
						high_jump_ratio[i] = -1.0; 
					}
			
					// reset adaptive counts and scale
					begin_jumps[i] = n_accepted[i]; 
					begin_draws[i] = n_draws; 
					end_draws[i] = n_draws+periods[i]; 
					scale[i] = new_scale; 
				}

				// not done?
				if (periods[i] <= max_period)
					done = false; 
			}
			check += period; 
		}
	}
	
	blocks.resize(k); 
	for (unsigned int i=0; i<k; i++)
		blocks[i] = B[i]*scale[i]; 
}

bool CMetropolis:: BlockRandomWalkMetropolis(double &log_posterior_y, TDenseVector &y, const TDenseVector &initial_v)
// Given
//	initial_v:	initial value, vector
// Results:
// 	y:	new vector according to random-walk metropolis
// 	log_posterior_y:	log_posterior(y)
// Returns:
// 	true:	if y is a new sample (different than x)
// 	false:	if y is the same as x
// Involving:
// 	blocks:	k matrices, each of which is n*b(i) consisting of the direction vectors of the i-th block
// 	model:	contains target distribution	
{
	size_t k=blocks.size(); 	// number of blocks	
	vector<size_t> b; 	// size of each block
	for (unsigned int i=0; i<k; i++)
		b[i] = blocks[i].cols; 

	TDenseVector x, y; 
	x.CopyContent(initial_v); 
	double log_previous = model->LogPosterior(x), log_current; 
	bool if_new_sample = false; 
	for (unsigned int i=0; i<k; i++)
	{
		y = x+blocks[i]*RandomNormalVector(b[i]); 
		log_current = model->LogPosterior(y); 
		if (log_current - log_previous >= log(dw_uniform_rnd()) )
		{
			x = y; 
			log_previous = log_current; 
			if (!if_new_sample)
				if_new_sample = true; 
		}
	}
	
	log_posterior_y = log_previous; 
	return if_new_sample; 
}
