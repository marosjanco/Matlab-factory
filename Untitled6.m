% Matlab code for plotting:

sol = csvread('3DGauss_Smooth2.csv');
NNN = length(sol);
B = round(NNN^(1/3));

ll = linspace(0,1,B+2);
[x,y,z] = meshgrid(ll,ll,ll);
val = reshape(sol,[B*B,B]);
V = zeros(size(x));

for i = 2 : B+1
    fix_z_soln = zeros(B+2,B+2);
    fix_z_soln(2:B+1,2:B+1) = reshape(val(:,i-1),[B,B]);
    V(:,:,i) = fix_z_soln;
end


figure
grid on
contourslice(x,y,z,V,[0.39],[0.39],[0.39])
xlim([0 1]); ylim([0 1]); zlim([0 1])
colorbar; axis equal
xlabel('x');ylabel('y');zlabel('z');
view(3)

figure
grid on
slice(x,y,z,V,[0.39],[0.39],[0.39])
colorbar;  axis equal
xlabel('x');ylabel('y');zlabel('z');
view(3)