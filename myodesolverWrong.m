function [xvec yvec] = myodesolverWrong( f, xrange, yinitial )
% calculate step size
xstep = (xrange(2) - xrange(1))/40;

% initialize values of x and y
x = xrange(1);
y = yinitial;
% initialize x-vector as [xrange(1) 0 0 ... 0]
xvec = [x; zeros(40,1)];
% initialize x-vector as [yinitial 0 0 ... 0]
yvec = [y; zeros(40,1)];
fprintf('xvec(1) = %f\n',xvec(1));
fprintf('yvec(1) = %f\n',yvec(1));
disp(' ');
% loop 41 times...
for iteration = 1:41
% calculate gradient
gradient = f(x,y);
% update value of y
y = y + gradient*xstep;
% update value of x
x = x + xstep;
% drop current x-value into x-vector
xvec(iteration) = x;
% drop current y-value into y-vector
yvec(iteration) = y;
fprintf('xvec(%d) = %f\n',iteration,xvec(iteration));
fprintf('yvec(%d) = %f\n',iteration,yvec(iteration));
disp(' ');
end
end

