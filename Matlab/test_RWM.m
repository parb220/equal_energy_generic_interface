% Code to test random walk metropolis code

rng('default');

n=1;
T=10000;
burnin=T/10;
thin=10;
x=zeros(n,1);
scale=20*eye(n);

mu=3*randn(n,1);
X=rand(n,n);
variance=X*X';
inverse_variance=inv(variance);

log_kernel=@(x)log_normal_kernel(x,mu,inverse_variance);

% test RandomWalkMetropolisSingle_dw()
% disp('Testing RandomWalkMetropolisSingle_dw()');
% X=zeros(n,T);
% log_current=log_kernel(x);
% for i=1:T
%     [x,log_current]=RandomWalkMetropolisSingle_dw(x, log_kernel, log_current, scale);
%     X(:,i)=x;
% end

% test MetropolisHastingsSingle_dw()
% disp('Testing MetropolisHastingsSingle_dw()');
% draw_proposal=@(x)RW_proposal(x,scale);
% X=zeros(n,T);
% [x log_proposal_kernel]=draw_proposal(x);
% log_ratio=log_kernel(x) - log_proposal_kernel;
% for i=1:T
%     [x,log_ratio]=MetropolisHastingsSingle_dw(x, log_kernel, draw_proposal, log_ratio);
%     X(:,i)=x;
% end

% test RandomWalkMetropolis_dw()
disp('Testing RandomWalkMetropolis_dw()');
[X, log_kernel_values acceptance]=RandomWalkMetropolis_dw(x,log_kernel,scale,T,burnin,thin);


% one dimensional test
z=X(1,:)';
I=eye(n);

disp('Actual mean and standard deviation');
mu(1)
sigma=sqrt(I(1,:)*variance*I(:,1))

disp('Mean and Standard Deviation');
sum(z')/T
sqrt((sum(z.*z) - (sum(z)^2)/T)/T)

acceptance

PlotCumulativeNormal_dw(z,mu(1),sigma);