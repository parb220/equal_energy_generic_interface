function blocks = BlockAdaptive(y,B,log_kernel,mid,period,max_period,verbose)
% y: n x 1 column vector - initial value
% B: k x 1 cell of n x b(i) matrices
% log_kernel: n-dimensional density kernel
% mid: target acceptance ratio
% period: initial number of draws to make before recomputing scale.  
%         increases by a factor of 2 when target ratio has been reached.
% max_period: final number of draws to make before recomputing scale.  
%             should be a power of 2 times period.
%
% verbose: 0 suppresses output
%
% Usually, b(1) + ... + b(k) = n and the matrix [B{1} ... B{k}] has 
% orthogonal columns, though this is not checked.
%
% Makes blockwise random-walk Metropolis draws, adjusting the block scale
% until the target acceptance rate is hit.  For each draw, the block i
% jumping kernel is as follows:
%
%   1) draw z from the b(i)-dimensional standard normal distribution
%   2) let x = scale(i)*B{i}*z
%   3) if y is the current point, then y + x is the proposal point
%

k=size(B,1);
n=size(y,1);

n_draws=0;
n_accepted=zeros(k,1);

previous_ratio=zeros(k,1);
scale=ones(k,1);
best_scale=ones(k,1);
low_scale=-ones(k,1);
low_jump_ratio=-ones(k,1);
high_scale=-ones(k,1);
high_jump_ratio=-ones(k,1);
begin_draws=zeros(k,1);
begin_jumps=zeros(k,1);
periods=period*ones(k,1);
end_draws=periods;

log_mid=log(mid);
lower_bound=exp(log_mid/0.2);
upper_bound=exp(log_mid/5.0);

b=zeros(k,1);
for i=1:k
    b(i)=size(B{i},2);
end

if verbose
    fprintf('Beginning adaptive burn-in -- maximum period = %d\n',max_period);
    fprintf('total ratio | current ratio | scale | period\n');
end

log_previous=log_kernel(y);
done=0;
check=period;
tic
while ~done
    % draw metropolis blocks
    for i=1:k
        x=y+scale(i)*B{i}*randn(b(i),1);
        log_current=log_kernel(x);

        if log_current - log_previous >= log(rand)
            y=x;
            log_previous=log_current;
            n_accepted(i)=n_accepted(i)+1;
        end
    end
    n_draws=n_draws+1;
    
    % rescale
    if n_draws == check
        if verbose
            fprintf('%d iterations completed\n',n_draws);
        end
        
        done=1;
        for i=1:k    
            % compute jump ratio
            previous_ratio(i) = (n_accepted(i) - begin_jumps(i))/(n_draws - begin_draws(i));
            
            % recompute scale?
            if (end_draws(i) <= n_draws) 
                
                % set new low or high bounds
                if previous_ratio(i) < mid
                    low_scale(i)=scale(i);
                    low_jump_ratio(i)=previous_ratio(i);
                else
                    high_scale(i)=scale(i);
                    high_jump_ratio(i)=previous_ratio(i);
                end
                
                % compute new scale and best scale
                 if low_jump_ratio(i) < 0.0
                     best_scale(i)=scale(i);
                     if previous_ratio(i) > upper_bound
                         new_scale=5.0*high_scale(i);
                     else
                         new_scale=(log_mid/log(previous_ratio(i)))*high_scale(i);
                     end
                 elseif high_jump_ratio(i) < 0.0
                     best_scale(i)=scale(i);
                     if previous_ratio(i) < lower_bound
                         new_scale=0.2*low_scale(i);
                     else
                         new_scale=(log_mid/log(previous_ratio(i)))*low_scale(i);
                     end
                 else
                     diff=high_jump_ratio(i) - low_jump_ratio(i);
                     if diff > 1.0e-6
                         new_scale=((mid - low_jump_ratio(i))*low_scale(i) + (high_jump_ratio(i) - mid)*high_scale(i))/diff;
                     else
                         new_scale=0.5*(low_scale(i) + high_scale(i));
                     end
                     best_scale(i)=new_scale; 
                     periods(i)=2*periods(i);
                     low_jump_ratio(i)=-1.0;
                     high_jump_ratio(i)=-1.0;
                 end           
                 
                 % reset adaptive counts and scale
                 begin_jumps(i)=n_accepted(i);
                 begin_draws(i)=n_draws;
                 end_draws(i)=n_draws + periods(i);
                 scale(i)=new_scale;
                 
                 % print output
                 if verbose
                     total_jump_ratio=n_accepted(i)/n_draws;
                     fprintf('block: %d  (%f %f %f %d)\n',i,total_jump_ratio,previous_ratio(i),scale(i),periods(i));
                 end                
            end
            
            % not done?
            if (periods(i) <= max_period)
                done=0;
            end
            
        end
        
        check=check+period;
    end
   
end
toc

blocks=cell(k,i);
for i=1:k
    blocks{i}=scale(i)*B{i};
end
