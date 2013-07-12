% four pass adaptive test

% set seed so parameters will be the same across runs
rng('default');

% set sizes
n=10;
ndraws=10000;
burnin=T/10;
thin=4;

% number Gaussians
m=4;

% random weights
w=rand(m,1);
w=w./sum(w);
log_w=log(w);

% mean, inverse variance, and log(abs(det(inverse variances)))
mu=randn(n,m);
iv=zeros(n,n,m);
log_det=zeros(m,1);

for i=1:m
    % get random orthogonal matrix Q
    [Q R]=qr(randn(n,n));
    
    % form diagonal
    d=0.1*(ones(n,1) + 0.1*(rand(n,1) - 0.5*ones(n,1)));
    
    % inverse variance
    iv(:,:,i)=Q*diag(1.0./d)*Q';
    
    % log(abs(det(inverse variance)))
    log_det(i)=-0.5*sum(log(d));
end

% set log kernel for target density
log_kernel=@(x)NormalMixture(x,log_w,mu,iv,log_det);

% set seed so runs will be different
rng('shuffle');

% run 4-pass adaptive
x=zeros(n,1);
period=100;
max_period=8*period;
block = FourPassAdaptive(x,log_kernel,period,max_period,verbose,ndraws,burnin,thin);

% simulate
[X, log_kernel_values acceptance]=BlockRandomWalkMetropolis_dw(x,log_kernel,block,T,burnin,thin);

acceptance

plot(log_kernel_values);
