%Cone
r=linspace(0,1,20);
theta0 =linspace(0,2*pi,40);
[r,theta0]=meshgrid(r,theta0);
x=r.*cos(theta0);
y=r.*sin(theta0);
z=r;
mesh(x,y,z)

%Circle
 radius = linspace(0,1,10); % For ten rings
 theta = (pi/180)*[0:15:360]; % For eight angles
 [R,T] = meshgrid(radius,theta); % Make radius/theta grid
 X = R.*cos(T); % Convert grid to cartesian coordintes
 Y = R.*sin(T);
 Z = 1+0.*X; %(corresponding values of Z in terms of X and Y)
 
 box on;
 grid on;
hold on;
mesh(X,Y,Z) %Circle
mesh(x,y,z)  %CONE

hold off;