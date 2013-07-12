% test MarginalDataDensity_MHM

% common values across runs
rng('default');

% sizes
n=2;
ndraws=10000;
burnin=ndraws/10;
thin=10;

% %=== normal posterior =====================================================
% mu=2*randn(n,1);
% X=randn(n,n);
% variance=X*X';
% inverse_variance=inv(variance);
% log_kernel=@(x)log_normal_kernel(x,mu,inverse_variance);
% 
% actual_log_mdd = -(-0.5*n*log(2*pi) + 0.5*log(abs(det(inverse_variance))));
% %==========================================================================

%=== mixture of normals posterior =========================================
% number Gaussians
m=4;

% random weights
w=rand(m,1);
w=w./sum(w);
log_w=log(w)

% mean, inverse variance, and log(abs(det(inverse variances)))
mu=randn(n,m);
iv=zeros(n,n,m);
log_det=zeros(m,1);

for i=1:m
    % get random orthogonal matrix Q
    [Q R]=qr(randn(n,n));
    
    % form diagonal
    d=1.0*(ones(n,1) + 0.5*(rand(n,1) - 0.5*ones(n,1)));
    
    % inverse variance
    iv(:,:,i)=Q*diag(1.0./d)*Q';
    
    % log(abs(det(inverse variance)))
    log_det(i)=-0.5*sum(log(d));
end

% set log kernel for target density
log_kernel=@(x)NormalMixture(x,log_w,mu,iv,log_det);

actual_log_mdd = 0.5*n*log(2*pi);
%==========================================================================

% different values across runs
rng('shuffle');

% setup single block adaptive procedure
period=200;
max_period=8*period;
verbose=1;
blocks = FourPassAdaptive(zeros(n,1),log_kernel,period,max_period,verbose,ndraws,burnin,thin);

%==========================================================================
%==========================================================================
%==========================================================================
% Setu USE_MODE to one to center at mode and zero to center at mean
USE_MODE=1;
    
% draw from posterior to setup proposal
[X, log_kernel_values acceptance]=BlockRandomWalkMetropolis_dw(zeros(n,1),log_kernel,blocks,ndraws,burnin,thin);

% use mode of the draws for the proposal
if USE_MODE == 1
    [s, i]=max(log_kernel_values');
    center=X(:,i);
else
    center=mean(X,2);
end

% variance
variance=cov(X');
inverse_variance_proposal=inv(variance);

% setup truncated log normal proposal
truncation=4.0;
constant=LogTruncatedNormalDensity(zeros(n,1),center,inverse_variance_proposal,truncation);
log_proposal_density=@(x)LogTruncatedNormalDensity(x,center,inverse_variance_proposal,truncation,constant);

% new draws of posterior
[X, log_kernel_values acceptance]=BlockRandomWalkMetropolis_dw(X(:,ndraws),log_kernel,blocks,ndraws,burnin,thin);

% compute log marginal data density (modified harmonic mean method
log_mdd_mhm = LogMarginalDataDensity_MHM(X, log_kernel_values, log_proposal_density)

actual_log_mdd

% compute log marginal data density (modified harmonic mean method)
scale=chol(variance,'lower');
draw_proposal=@()DrawNormal(center,scale);
log_proposal_density=@(x)log_normal_density(x,center,inverse_variance_proposal);
n_proposal_draws=10000;
log_mdd_bridge = LogMarginalDataDensity_Bridge(X, log_kernel_values, log_kernel, draw_proposal, log_proposal_density, n_proposal_draws)

% compute log marginal data density (modified harmonic mean method)
scale=chol(variance,'lower');
draw_proposal=@()DrawNormal(center,scale);
log_proposal_density=@(x)log_normal_density(x,center,inverse_variance_proposal);
n_proposal_draws=10000;
log_mdd_wz = LogMarginalDataDensity_WZ(X, log_kernel_values, log_kernel, draw_proposal, log_proposal_density, n_proposal_draws)
