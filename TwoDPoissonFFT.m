TwoDsol = csvread('2DSmoothFastSN.csv');
x = linspace(0,1,length(TwoDsol));
y = linspace(0,1,length(TwoDsol));
[X,Y] = meshgrid(x,y);
contour(X,Y,TwoDsol);
colorbar;
xlabel('X')
ylabel('Y')
axis equal