function  [Q2bc, interfW] = setBoundary2p(q1,N1,M1,Q2bc,N2,M2)

        global dx

        sol1 = reshape(q1,N1,M1);
        interfW = (sol1(48,:)+sol1(49,:))/2;
        Q2bc           = reshape(Q2bc,N2,M2);
        Q2bc(1,33:end) = interfW*2/dx/dx;
        Q2bc           = Q2bc(:);
end