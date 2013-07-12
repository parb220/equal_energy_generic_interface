function x = AddLogs(y,z)
% y = log(u)
% z = log(v)
% 
% x = log(u + v)

if y > z
    x = y + log(exp(z - y) + 1.0);
else
    x = z + log(1.0 + exp(y - z));
end