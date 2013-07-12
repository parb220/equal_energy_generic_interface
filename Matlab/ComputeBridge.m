function log_mdd = ComputeBridge(proposal,posterior)
%
% proposal:
%   n2 x 2 matrix.  First column is log proposal density and second column
%   is log posterior kernel of proposal draws
%
% posterior:
%   n1 x 2 matrix.  First column is log proposal density and second column
%   is log posterior kernel of posterior draws
%

MAX_C = 1.0e50;
TOL = 1.0e-7;

% Bracket the zero
diff=BridgeDifference(proposal,posterior,0.0);
if diff < 0.0
    max_c=0.0;
    min_c=-1.0;
    diff=BridgeDifference(proposal,posterior,min_c);
    while (diff <= 0) && (min_c > -MAX_C)
        max_c=min_c;
        min_c=10*min_c;
        diff=BridgeDifference(proposal,posterior,min_c);
    end
else
    min_c=0.0;
    max_c=1.0;
    diff=BridgeDifference(proposal,posterior,max_c);
    while (diff >= 0) && (max_c < MAX_C)
        min_c=max_c;
        max_c=10*max_c;
        diff=BridgeDifference(proposal,posterior,max_c);     
    end
end

% Divide and conququer
if max_c >= MAX_C
    log_mdd=max_c;
elseif min_c <= -MAX_C
    log_mdd=min_c;
else
    mid_c=0.5*(min_c + max_c);
    diff=BridgeDifference(proposal,posterior,mid_c);
    i=0;
    while (i < 50) && (abs(diff) >= TOL)
        if diff > 0
            min_c=mid_c;
        else
            max_c=mid_c;
        end
        mid_c=0.5*(min_c + max_c);
        diff=BridgeDifference(proposal,posterior,mid_c);
    end
    log_mdd=mid_c;
end