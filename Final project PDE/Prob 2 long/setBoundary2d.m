function  [Q2bc, interfW] = setBoundary2d(q1,N1,M1,Q2bc,N2,M2,q2Nfirst)

        global Re dy dx dt

        sol1 = reshape(q1,N1,M1);
        interfW = (sol1(48,:)+sol1(49,:))/2;
        Q2bc           = reshape(Q2bc,N2,M2);
        Q2bc(1,33:end) = interfW*2*dt/Re/dx/dx;
        Q2bc(1,end)    = Q2bc(1,end) + q2Nfirst*2*dt/Re/dy/dy;
        Q2bc           = Q2bc(:);
        
end