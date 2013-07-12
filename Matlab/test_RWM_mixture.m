% Code to test random walk metropolis code

rng('default');

n=2;
T=10000;
burnin=T/10;
thin=10;

scale=0.09*eye(n);

m=2;
w=rand(m,1);
w=w./sum(w);
log_w=log(w);
mu=zeros(n,m);
iv=zeros(n,n,m);
log_det=zeros(m,1);

for i=1:m
    mu(i,i)=i;
    
    [Q R]=qr(randn(n,n));
    d=.01*(ones(n,1) + 0.0*(rand(n,1) - 0.5*ones(n,1)));
    iv(:,:,i)=Q*diag(1.0./d)*Q';
    log_det(i)=-0.5*sum(log(d));
end
    
    
% test RandomWalkMetropolis_dw()
rng('shuffle');
disp('Testing RandomWalkMetropolis_dw()');
log_kernel=@(x)NormalMixture(x,log_w,mu,iv,log_det);
x=zeros(n,1);
[X, log_kernel_values acceptance]=RandomWalkMetropolis_dw(x,log_kernel,scale,T,burnin,thin);

acceptance

mu

scatter(X(1,:),X(2,:));


