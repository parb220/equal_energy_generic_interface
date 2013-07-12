% test block adaptive

% set seed so parameters will be the same across runs
rng('default');

% set sizes
n=2;
T=10000;
burnin=T/10;
thin=10;

% number Gaussians
m=1;

% random weights
w=rand(m,1);
w=w./sum(w);
log_w=log(w);

% mean, inverse variance, and log(abs(det(inverse variances)))
mu=zeros(n,m);
iv=zeros(n,n,m);
log_det=zeros(m,1);

for i=1:m
    % mean - requires m <= n
    mu(i,i)=i;
    
    % get random orthogonal matrix Q
    [Q R]=qr(randn(n,n));
    
    % form diagonal
    d=0.1*(ones(n,1) + 0.5*(rand(n,1) - 0.5*ones(n,1)));
    
    % inverse variance
    iv(:,:,i)=Q*diag(1.0./d)*Q';
    
    % log(abs(det(inverse variance)))
    log_det(i)=-0.5*sum(log(d));
end

% set log kernel for target density
log_kernel=@(x)NormalMixture(x,log_w,mu,iv,log_det);

% set seed so runs will be different
rng('shuffle');

% setup block adaptive
n_blocks=2;
B=cell(n_blocks,1);
y=zeros(n,1);
period=1000;
max_period=8*period;
verbose=1;

if n_blocks == 1
    % a single block
    B{i}=eye(n);
elseif n_blocks == 2
    % two blocks of almost equal size
    I=eye(n);
    m=floor(n/2);
    B{1}=I(:,1:m);
    B{2}=I(:,m+1:end);
end

blocks = BlockAdaptive(y,B,log_kernel,0.25,period,max_period,verbose);
    
% test RandomWalkMetropolis_dw()
disp('Testing RandomWalkMetropolis_dw()');

[X, log_kernel_values acceptance]=BlockRandomWalkMetropolis_dw(y,log_kernel,blocks,T,burnin,thin);

acceptance

% one dimensional test - assumes number of Gaussians (m) is one.
if m == 1
    z=X(1,:)';
    I=eye(n);
    
    disp('Actual mean and standard deviation');
    mu(1)
    sigma=sqrt(I(1,:)*inv(iv(:,:,1))*I(:,1))
    
    disp('Mean and Standard Deviation');
    sum(z')/T
    sqrt((sum(z.*z) - (sum(z)^2)/T)/T)
    
    PlotCumulativeNormal_dw(z,mu(1),sigma);
end





