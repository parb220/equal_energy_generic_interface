function [x log_current] = MetropolisHastingsSingle_dw(y, log_target_kernel, draw_proposal, log_previous)
% y                   - current draw
% log_target_kernel   - function to compute the posterior kernel, takes a 
%                       single n-dimensional vector argument.
% proposal_draw       - function make draw from proposal distribution, 
%                       takes a single n-dimensional vector argument, which
%                       is the current draw.  Returns draw and the log of 
%                       proposal kernel.
% log_previous        - log_target_kernel(y) - log_proposal_kernel(y)
%
% Returns a n x 1 Metropolis-Hastings draw from the target distribution 
% and the log ratio of the kernels.

[x log_proposal_kernel]=draw_proposal(y);
log_current=log_target_kernel(x) - log_proposal_kernel;

if log_current - log_previous < log(rand)
    x=y;
    log_current=log_previous;
end
