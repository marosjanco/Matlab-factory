[X,Y] = meshgrid(-1.5:.05:1.5);
Z = (X-Y).*(X.^2+Y.^2-1);

figure
mesh(X,Y,Z)
axis([-1 1 -1 1 -1 1])
view(-45,10)

figure
contour(X,Y,Z,35)