% previous value + gradient * (x - x0)


% 
% gamma = 1.4;
% pinitial = 2000;
% pfinal = 100;
% pstep = -100;
% m = 40;
% 
% 
% for p = pinitial:pstep:pfinal-pstep
% 
% gradient = 1/gamma * m/p;
% m = m + gradient*pstep;
% 
% end
% m


function [xvec, yvec] = myodesolver( f, xrange, yinitial )

xstep = (xrange(2) - xrange(1))/40;

x = xrange(1);
y = yinitial;

xvec = [x; zeros(40,1)];
yvec = [y; zeros(40,1)];

for iteration = 1:40

gradient = f(x,y);
y = y + gradient*xstep;
x = x + xstep;

xvec(iteration+1) = x;
yvec(iteration+1) = y;
end
end

% function [xvec yvec] = myodesolver( f, xrange, yinitial )
% 
% [pvals mvals] = myodesolver(@(p, m) 1/1.4*m/p, [2000 100], 40)
% plot(pvals, mvals);
% xlabel('P/Pa');
% ylabel('m/kg');
% title('Plot of mass m against pressure P for a pressure release');