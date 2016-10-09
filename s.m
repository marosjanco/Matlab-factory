function [S] = s(n) 
S = @(x)0;
for m = 0:n
    N = @(x)(1./(n+1));
    F_0 = fs(m);
    F = @(x) N(x).*F_0(x);
    S = @(x)S(x)+F(x);
end
end
% f = @(x) x ;
% g = @(x) x^2 ;
% h = @(x) f(x) + g(x) 