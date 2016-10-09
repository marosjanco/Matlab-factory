sol = csvread('3DGauss_Smooth.csv');
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
grid on;
    contourslice(x,y,z,V,[0.4],[0.4],[0.4]); 
    slice(x,y,z,V,[0.4],[0.4],[0.4]);
    colorbar;
    view(3);
    hold off;