function [X log_kernel acceptance] = RandomWalkMetropolis_dw(x, log_target_kernel, scale, ndraws, burnin, thin)
% x                 - initial value (n x 1) vector
% log_target_kernel - function to compute the posterior kernel, takes a 
%                     single n-dimensional vector argument
% scale             - sqare root of variance-covariance matrix for normal 
%                     jumping proposal
% ndraws            - number of draws to retain
% burnin            - initial draws to discard
% thin              - total draws made burnin + ndraws * thin
% 
% Return n x ndraws matrix of random walk Metropolis draws from the target 
% distribution.

n=size(x,1);
X=zeros(n,ndraws);
log_kernel=zeros(1,ndraws);

log_previous=log_target_kernel(x);
a=0;

for i=1:burnin
    
    z=x+scale*randn(n,1);
    log_current=log_target_kernel(z);
    
    if log_current - log_previous >= log(rand)
        log_previous=log_current;
        x=z;
        a=a+1;
    end
    
end

for i=1:ndraws
    
    for j=1:thin
        
        z=x+scale*randn(n,1);
        log_current=log_target_kernel(z);
        
        if log_current - log_previous >= log(rand)
            log_previous=log_current;
            x=z;
            a=a+1;
        end
    end
    
    X(:,i)=x;
    log_kernel(i)=log_previous;
    
end

acceptance = a / (burnin + ndraws*thin);