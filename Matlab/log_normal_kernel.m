function ln_p = log_normal_kernel(x, mu, inverse_variance)

ln_p = -0.5 * (x - mu)' * inverse_variance * (x - mu);