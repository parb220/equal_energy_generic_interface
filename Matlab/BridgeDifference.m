function difference = BridgeDifference(proposal, posterior, logr)
%
% proposal:
%   n2 x 2 matrix.  First column is log proposal density and second column
%   is log posterior kernel of proposal draws
%
% posterior:
%   n1 x 2 matrix.  First column is log proposal density and second column
%   is log posterior kernel of posterior draws
%
% logr
%   test value of log mdd
%

MINUS_INFINITY = -1.0e300;
n1=size(posterior,1);
n2=size(proposal,1);

r=n2/n1;
sum2=MINUS_INFINITY;
for i=1:n2
    x=proposal(i,1) + logr - proposal(i,2);
    if x < 0
      sum2=AddLogs(sum2,-log(1+r*exp(x)));
    else
      sum2=AddLogs(sum2,-x-log(exp(-x)+r));
    end
end

r=n1/n2;
sum1=MINUS_INFINITY;
for i=1:n1
    x=posterior(i,2) - posterior(i,1) - logr;
    if x < 0
      sum1=AddLogs(sum1,-log(1+r*exp(x)));
    else
      sum1=AddLogs(sum1,-x-log(exp(-x)+r));
    end
end
    
difference=sum2 - sum1;