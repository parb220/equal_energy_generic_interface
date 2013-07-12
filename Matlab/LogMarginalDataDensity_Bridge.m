function log_mdd = LogMarginalDataDensity_Bridge(posterior_draws, log_kernel_values, log_kernel, draw_proposal, log_proposal_density, n_proposal_draws)
% draws:
%   n x T matrix containing the posterior draws in the columns. 
%
% log_kernel_values: 
%   1 x T column vector containing the log posterior kernal values of the 
%   draws.
%
% log_proposal_density: 
%   Function that computes the proposal density.  Takes a single
%   n-dimensional column vector argument. 
%
% Returns the marginal data density.
%

n_posterior_draws=size(posterior_draws,2);
posterior=zeros(n_posterior_draws,2);
for i=1:n_posterior_draws
    posterior(i,1)=log_proposal_density(posterior_draws(:,i));
    posterior(i,2)=log_kernel_values(i);
end

proposal=zeros(n_proposal_draws,2);
for i=1:n_proposal_draws
    x=draw_proposal();
    proposal(i,1)=log_proposal_density(x);
    proposal(i,2)=log_kernel(x);
end

log_mdd=ComputeBridge(proposal,posterior);