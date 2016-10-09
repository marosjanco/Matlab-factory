function cobweb( f, x0, nIterations, xmin, xmax )


xvec = zeros(2,2*nIterations+1);

xvec(:,1) = [x0;x0];


oldx = x0;

for iteration = 1:nIterations

newx = f(oldx);
xvec(:,2*iteration) = [oldx;newx];
xvec(:,2*iteration+1) = [newx;newx];
oldx = newx;
end

plot(xvec(1,:),xvec(2,:),'r');


hold on;
xvals = linspace(xmin,xmax,20);
yvals = f(xvals);
plot(xvals,xvals,'k');
plot(xvals,yvals,'b');
hold off;

axis equal;
axis([xmin, xmax, xmin, xmax]);
end