% normal importance sampling

n=1;
N=100000;

% set generator to have repeatable variances
rng('default');

% set up variance 
[Q R]=qr(randn(n,n));
d=ones(n,1) + 1.5*(rand(n,1) - 0.5*ones(n,1));
variance=Q*diag(d)*Q';
inverse_variance=Q*diag(1.0./d)*Q';
sqrt_det=sqrt(prod(d));

% set up mean
mu=3.0*sqrt(n)*ones(n,1);

rng('shuffle');

% importance sampler
y=zeros(1+n,N);
for i=1:N
    
    y(2:n+1,i) = randn(n,1);
     
    y(1,i) = 1.0/sqrt_det * exp(-0.5*(y(2:n+1,i) - mu)'*inverse_variance*(y(2:n+1,i) - mu)) / exp(-0.5*y(2:n+1,i)'*y(2:n+1,i));
    
end

% tests
disp('Max weight');
max(y(1,:))

disp('Mean weight');
mean(y(1,:))

disp('Mean');
m=zeros(n,1);
sum=0;
for i=1:N
    m=m + y(1,i)*y(2:n+1,i);
    sum=sum + y(1,i);
end
[m/sum mu]

plot(y(1,:));
