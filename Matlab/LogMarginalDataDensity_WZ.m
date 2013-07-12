function log_mdd = LogMarginalDataDensity_WZ(posterior_draws, log_kernel_values, log_kernel, draw_proposal, log_proposal_density, n_proposal_draws)
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

S1=sort(proposal(:,2),'descend');
if size(S1,1) > 30
    a1=S1(30);
else
    a1=S1(end);
end
S2=sort(posterior(:,2),'descend');
a2=S2(floor(size(S2,1)/10));

if a1 >= a2
    L1=a2;
else
    L1=a1;
end
    
[log_mdd, in_proposal, in_posterior, I]=ComputeMDD_WZ(proposal,posterior,L1,1.0e300);

fprintf('%d out of %d proposal draws in region\n',in_proposal,size(proposal,1));
fprintf('%d out of %d posterior draws in region\n',in_posterior,size(posterior,1));

