function interpolation(a,b,n)
% function sheet7q1(a,b,n)
%
% Interpolate f(x) = 1 / (1 + x^2) on [a,b] with equally spaced points and
% Chebyshev points. Plot f, p_n, p_n^c and errors f-p_n, f-p_n^c.

xc = ones(1,n+1);
x = (a:(b-a)/n:b);
for i = 0:n
xc(i+1) = 0.5 * ( (b-a) * cos( (2*i+1) / (n+1) * pi / 2 ) + a + b );
end
f = 1 ./ ( 1 + x.*x); p = polyfit(x, f, n);
fc = 1 ./ ( 1 + xc.*xc); pc = polyfit(xc, fc, n);
xx = (a:(b-a)/1000:b); f = 1 ./ ( 1 + xx.*xx);
close all;
plot(xx,f,'k',xx,polyval(p,xx),'r',xx,polyval(pc,xx),'b');
title('f, p_n, p_n^c');
figure;
subplot(2,1,1); plot(xx,f-polyval(p,xx),'r',xx,f-polyval(pc,xx),'b');
title('f - p_n, f - p_n^c');
pro = ones(size(xx)); proc = pro; cbound = 2;
for i = 0:n
pro = pro .* (xx - x(i+1)) ./ (i+1);
proc = proc .* (xx - xc(i+1)) ./ (i+1);
cbound = cbound * (b-a) / 2 / (i+1) / 2
end
subplot(2,1,2); plot(xx,pro,'r',xx,proc,'b',xx,cbound,'k',xx,-cbound,'k');
title('1/(n+1)! \Pi_{i=0}^n (x-x_i)');
