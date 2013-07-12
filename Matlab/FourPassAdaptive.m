function block = FourPassAdaptive(x,log_kernel,period,max_period,verbose,ndraws,burnin,thin)
% a 4-pass adaptive technique
% first pass
%   identity variance matrix / n blocks
% second pass
%   diagonal variance matrix / 1 block
% third pass
%   simulate from second pass block structure
%   variance matrix from simulation / n blocks
% forth pass
%   variance matrix from simulation with scaled columns / 1 block

n=size(x,1);

% first pass - n blocks
n_blocks=n;
B=cell(n_blocks,1);
I=eye(n);
for i=1:n_blocks
    B{i}=I(:,i);
end
blocks = BlockAdaptive(x,B,log_kernel,0.25,period,max_period,verbose);

% second pass - 1 block
n_blocks=1;
B=cell(n_blocks,1);
B{1}=zeros(n,n);
for i=1:n
    B{1}(:,i)=blocks{i};
end
blocks = BlockAdaptive(x,B,log_kernel,0.25,period,max_period,verbose);

% simulate
[X log_kernel_values acceptance] = BlockRandomWalkMetropolis_dw(x,log_kernel,B,ndraws,burnin,thin);

% compute variance covariance
variance=zeros(n,n);
for i=1:ndraws
    variance=variance + X(:,i)*X(:,i)';
end
variance=0.5*(variance + variance')/ndraws;
[U D V]=svd(variance);
D=sqrt(D);
U=U*D;

% third pass - n blocks
n_blocks=n;
B=cell(n_blocks,1);
for i=1:n_blocks
    B{i}=U(:,i);
end
blocks = BlockAdaptive(x,B,log_kernel,0.25,period,max_period,verbose);

% forth pass - 1 block
n_blocks=1;
B=cell(n_blocks,1);
B{1}=zeros(n,n);
for i=1:n
    B{1}(:,i)=blocks{i};
end
block = BlockAdaptive(x,B,log_kernel,0.25,period,max_period,verbose);


