function [F] = fs(in) 
N = ceil(in/2);
if N==0 
    F = @(x)0;
else
    
F = @(x)0*x;
for m = 1:N
    K = @(x)(4./((2*m-1).*pi))*sin((2*m-1).*pi.*x./2);
    F = @(x)(F(x)+K(x));
end
end
end
