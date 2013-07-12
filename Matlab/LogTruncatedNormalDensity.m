function log_density = LogTruncatedNormalDensity(x, mu, inverse_variance, truncation, constant)
% x = n x 1 column vector at which to evaluate density
% mu = n x 1 column vector which is the mean of underlying normal
% inverse_variance = inverse of the variance of underlying normal
% truncation = truncation level in variance units
% constant = 0.5*log(2*pi) + 0.5*log(abs(det(inverse_varinace))) -
%                                                 cdf('chi2',truncation,n)
%
% Returns log of the truncated normal density.  If constant is not passed,
% then it is computed and returned in log_density.
%

n=size(x,1);

if nargin < 5
    log_density = -0.918938533204673*n + 0.5*log(abs(det(inverse_variance))) - log(cdf('chi2',truncation,n));
else
    v=(x - mu)'*inverse_variance*(x - mu);
    
    if v > truncation
        log_density=-1e300;
    else
        log_density = constant - 0.5*v;
    end
end

