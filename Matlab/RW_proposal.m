function [x log_kernel]=RW_proposal(y, scale)
% normal draw centered at y

z=randn(size(y));
log_kernel = -0.5*dot(z,z);
x=y + scale*z;
