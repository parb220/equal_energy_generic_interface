function x = AddScaledLogs(a,y,b,z)
% y = log(u)
% z = log(v)
% 
% x = log(a*u + b*v)

if y > z
    x = y + log(b*exp(z - y) + a);
else
    x = z + log(b + a*exp(y - z));
end