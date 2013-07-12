% One dimensional normal draws


N=50000;
N2=2*N;
x=zeros(N2,1);
y=zeros(N2,1);
z=zeros(N2,1);
mu=3.6;
sigma=2.4;

%==========================================================================
tic
for i=1:N2
    x(i)=mu + sigma*randn;
end
toc

disp('Mean and Standard Deviation');
sum(x)/N2
sqrt((sum(x.*x) - (sum(x)^2/N2))/N2)

%==========================================================================
tic
for i=1:N2
    y(i)=icdf('Normal',rand,mu,sigma);
end
toc

disp('Mean and Standard Deviation');
sum(y)/N2
sqrt((sum(y.*y) - (sum(y)^2)/N2)/N2)

%==========================================================================
tic
c=2.0*pi;
for i=1:N
    a=sqrt(-2.0*log(rand));
    u=c*rand;
    z(2*i-1)=mu + sigma*a*cos(u);
    z(2*i)=mu + sigma*a*sin(u);
end
toc

disp('Mean and Standard Deviation');
sum(z)/N2
sqrt((sum(z.*z) - (sum(z)^2)/N2)/N2)

%PlotCumulativeNormal_dw(z,mu,sigma);

%==========================================================================
tic
w=mu*ones(N2,1) + sigma*randn(N2,1);
toc

disp('Mean and Standard Deviation');
sum(w)/N2
sqrt((sum(w.*w) - (sum(w)^2/N2))/N2)

    