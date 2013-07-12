function ln_p = log_normal_density(x, mu, inverse_variance)

n=size(x,1);
ln_p = -0.918938533204673*n + 0.5*log(abs(det(inverse_variance))) - 0.5 * (x - mu)' * inverse_variance * (x - mu);