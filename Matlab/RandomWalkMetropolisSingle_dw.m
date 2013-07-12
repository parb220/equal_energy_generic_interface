function [x log_current] = RandomWalkMetropolisSingle_dw(y, log_target_kernel, log_previous, scale)
% y                 - current draw
% log_target_kernel - function to compute the posterior kernel, takes a 
%                     single n-dimensional vector argument.
% previous          - target_kernel(y)
% scale             - sqare root of variance-covariance matrix for normal 
%                     jumping proposal
%
% Returns a n x 1 random walk Metropolis draw from the target distribution 
% and its log kernel.

n=size(y,1);
x=y+scale*randn(n,1);
log_current=log_target_kernel(x);

if log_current - log_previous < log(rand)
    x=y;
    log_current=log_previous;
end
