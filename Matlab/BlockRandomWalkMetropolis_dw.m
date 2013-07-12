function [X log_kernel acceptance] = BlockRandomWalkMetropolis_dw(x, log_target_kernel, blocks, ndraws, burnin, thin)
% x                 - initial value (n x 1) vector
% log_target_kernel - function to compute the posterior kernel, takes a 
%                     single n-dimensional vector argument
% blocks            - k x 1 cell of n x b(i) matrices
% scale             - sqare root of variance-covariance matrix for normal 
%                     jumping proposal
% ndraws            - number of draws to retain
% burnin            - initial draws to discard
% thin              - total draws made burnin + ndraws * thin
% 
% Return n x ndraws matrix of random walk Metropolis draws from the target 
% distribution.
%
% Usually, b(1) + ... + b(k) = n and the matrix [B{1} ... B{k}] has 
% orthogonal columns, though this is not checked.
%
% Makes blockwise random-walk Metropolis draws.  For each draw, the block i
% jumping kernel is as follows:
%
%   1) draw z from the b(i)-dimensional standard normal distribution
%   2) let x = scale(i)*B{i}*z
%   3) if y is the current point, then y + x is the proposal point
%

n=size(x,1);
k=size(blocks,1);
X=zeros(n,ndraws);
log_kernel=zeros(1,ndraws);

log_previous=log_target_kernel(x);
a=zeros(k,1);
b=zeros(k,1);
for i=1:k
    b(i)=size(blocks{i},2);
end

fprintf('beginning %d burn-in draws',burnin);
tic
for m=1:burnin  
    % block random walk metropolis
    for i=1:k
        y=x+blocks{i}*randn(b(i),1);
        log_current=log_target_kernel(y);
        
        if log_current - log_previous >= log(rand)
            x=y;
            log_previous=log_current;
            a(i)=a(i)+1;
        end
    end
end
toc

total=ndraws*thin;
fprintf('beginning %d draws\n',total);
tic
for m=1:ndraws
    % thin
    for j=1:thin
        % block random walk metropolis
        for i=1:k
            y=x+blocks{i}*randn(b(i),1);
            log_current=log_target_kernel(y);
            
            if log_current - log_previous >= log(rand)
                x=y;
                log_previous=log_current;
                a(i)=a(i)+1;
            end
        end
    end
    
    X(:,m)=x;
    log_kernel(m)=log_previous;
    
end
toc

acceptance = a / (burnin + ndraws*thin);