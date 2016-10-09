 [x, y] = meshgrid(-1:0.05:2, -3:.05:2);
xdot = (1-y-x).*x; %Note the use of .* and .^
ydot = (x-2-y).*y;
quiver(x,y,xdot, ydot)
axis equal;
% 
% [x, y] = meshgrid(-.5:0.05:0.5, -.5:.05:.5);
% xdot = x - y; %Note the use of .* and .^
% ydot = x+3*y;
% quiver(x,y,xdot, ydot)
% axis equal;
% % % 
% % phase portraits in matlab
% 
% 
% [x,y]=meshgrid(-3:.3:3,-2:.3:2);
% dx= -x + 6*(y);
% dy = -x + (y);
% quiver(x,y,dx,dy)