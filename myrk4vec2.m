function [tvec, xvec] = myrk4vec2(f, nvector, x0)
format long
dim=length(x0);
 % calculate step size
tstep = range(nvector)/(length(nvector)-1);
% initialize values of t and x
t = nvector(1);
x = x0;
% initialize tvec as [nvector(1) 0 0 ... 0]
tvec = [t; zeros(length(nvector)-1,1)];
% initialize xvec as [x0 0 0 ... 0]
xvec = zeros(length(nvector),dim);
xvec(1,:)=x;
% loop...
for i = 1:length(nvector)-1
    % calculate k1,k2,k3,k4
    x=xvec(i,:);
    t=tvec(i);
    k1= tstep*(f(t,x))';
    k2= tstep*(f(t+tstep/2, x+k1/2))';
    k3= tstep*(f(t+tstep/2, x+k2/2))';
    k4= tstep*(f(t+tstep, x+k3))';
    % calculate gradient
    gradient = (k1+2*k2+2*k3+k4)/6;
    % update value of t and x
    x=(x+gradient)';
    t=t + tstep;
    % drop current t-value into t-vector
    tvec(i+1) = t;
    % drop current x-value into x-vector
    xvec(i+1,:) = x;
end
end
