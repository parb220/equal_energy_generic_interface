function log_kernel = NormalMixture(x, log_w, mu, iv, log_det)
% x       - n dimensional column vector
% log_w   - m dimensional column vector of log weights
% mu      - n x m matrix of means
% iv      - n x n x m matrix of inverse variances
% log_det - 0.5*log(abs(det(iv(:,:,i))))
%
% returns the log of the density of the mixture of normals
%

n=size(x,1);
m=size(log_w,1);

log_kernel=log_w(1) + log_det(1) - 0.5*(x-mu(:,1))'*iv(:,:,1)*(x-mu(:,1));

for i=2:m
    log_kernel=AddLogs(log_kernel,log_w(i) + log_det(i) - 0.5*(x-mu(:,i))'*iv(:,:,i)*(x-mu(:,i)));
end