function PlotCumulativeNormal_dw(y, mu, sigma)
% Plots cumulative distribution from x and actual cumulative distribution

T=size(y,1);
y=sort(y);

plot_points=100;

minimum=y(1);
maximum=y(end);
inc=(maximum - minimum)/(plot_points-1);

x=(minimum:inc:maximum)';

w=zeros(plot_points,1);
for i=1:plot_points
    w(i)=sum(y <= x(i))/T;
end
    
%z=[cdf('Normal',x,mu,sigma) w];

plot(x,cdf('Normal',x,mu,sigma),x,w);