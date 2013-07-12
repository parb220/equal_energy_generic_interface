function log_mdd = LogMarginalDataDensity_MHM(draws, log_kernel_values, log_proposal_density)
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

T=size(draws,2);

log_mdd=log_proposal_density(draws(:,1)) - log_kernel_values(1,1);

for i=2:T
    p=log_proposal_density(draws(:,i)) - log_kernel_values(1,i);
    log_mdd=AddLogs(log_mdd,p);
end

log_mdd=-log_mdd + log(T);