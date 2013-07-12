function [log_mdd, in_proposal, in_posterior, I] = ComputeMDD_WZ(proposal, posterior, L1, L2)
%    Let h(x) and f(x) be properly scaled probability density functions and 
%    let c be an unknown constant.  We use the following notation.
% 
%       x      - parameters
%       y      - data
%       f(x)   - posterior distribution = p(x|y)
%       h(x)   - proposal distribution
%       c      - marginal data distribution = p(y)
%       c*f(x) - likelihood*prior = p(y|x)*p(x)
% 
%    Assumes:
%      proposal  : N x 2 matrix with proposal(i,0) = ln(h(x(i))) and 
%                  proposal(i,1) = ln(c*f(x(i))) where x(i) is sampled from 
%                  h(x).
%      posterior : M x 2 matrix with posterior(i,0) = ln(h(x(i))) and 
%                  posterior(i,1) = ln(c*f(x(i))) where x(i) is sampled 
%                  from f(x).
%      L1        : cutoff value (c*f(x) > L1)
%      L2        : cutoff value for (h(x) < L2)
% 
%    Returns:
%      Estimate of c or MINUS_INFINITY if no proposal draws satisfied the
%      restriction given by the cutoff values of L1 and L2.
% 
%    Notes:
%      Let S be the set of all x such that c*f(x) > exp(L1) and 
%      h(x) < exp(L2). Then c = p(L1,L2)/I(L1,L2) where p(L1,L2) is the 
%      probability that x, sampled from h(x), is in S and I(L1,L2) is the 
%      integral with respect to x over the set S of the integrand
% 
%                             h(x)/(c*f(x)) * f(x) 
% 
%      p(L1,L2) can be approximated from the proposal draws by
% 
%                (number of draws in S) / (total number of draws)
% 
%      I(L1,L2) can be approximated from the posterior draws by summing 
% 
%                                     h(x)/(c*f(x)) 
%    
%      over all posterior draws with x in S and then dividing by the total 
%      number of posterior draws.       
%

MINUS_INFINITY=-1.0e300;
PLUS_INFINITY=1.0e300;

if (L1 <= MINUS_INFINITY)
    if (L2 >= PLUS_INFINITY)
        in_proposal=size(proposal,1);
        in_posterior=size(posterior,1);
        I=MINUS_INFINITY;
        for i=1:size(posterior,1)
            I=AddLogs(I,posterior(i,1) - posterior(i,2));
        end
    else
        in_proposal=0;
        for i=1:size(proposal,1)
            if proposal(i,1) <= L2
                in_proposal=in_proposal+1;
            end
        end
        
        I=MINUS_INFINITY;
        in_posteror=0;
        for i=1:size(posterior,1)
            if posterior(i,1) <= L2
                in_posterior=in_posteror+1;
                I=AddLogs(I,posterior(i,1) - posterior(i,2));
            end
        end
    end                         
else
    if L2 >= PLUS_INFINITY
        in_proposal=0;
        for i=1:size(proposal,1)
            if proposal(i,2) >= L1
                in_proposal=in_proposal+1;
            end
        end

        I=MINUS_INFINITY;
        in_posterior=0;
        for i=1:size(posterior,1)
            if posterior(i,2) >= L1
                in_posterior=in_posterior+1;
                I=AddLogs(I,posterior(i,1) - posterior(i,2));
            end
        end
    else
        in_proposal=0;
        for i=1:size(proposal)
            if (proposal(i,2) >= L1) && (proposal(i,1) <= L2)
                in_proposal=in_proposal+1;
            end
        end
        
        I=MINUS_INFINITY;
	    in_posterior=0;
        for i=1:size(posterior,1)
            if (posterior(i,2) >= L1) && (posterior(i,1) <= L2)
                in_posterior=in_posterior+1;
	            I=AddLogs(I,posterior(i,1) - posterior(i,2));
            end
        end
    end
end

if in_posterior > 0
    I=I-log(size(posterior,1));
end

if in_proposal > 0
    log_mdd=log(in_proposal/size(proposal,1)) - I;
else
    log_mdd=MINUS_INFINITY;
end
