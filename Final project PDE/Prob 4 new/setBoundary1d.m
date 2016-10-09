function [Q1bc, interfS, interfE] = setBoundary1d(q2,N2,M2,Q1bc,N1,M1,q1Nlast)
        
        global Re dy dx dt
        sol2    = reshape(q2,N2,M2);
        
        interfS = (sol2(1:16,32)+sol2(1:16,33))/2;
        interfE = (sol2(16,33:end)+sol2(17,33:end))/2;
      
        Q1bc             = reshape(Q1bc,N1,M1);
        Q1bc(end,:)      = interfE*2*dt/Re/dx/dx;
        Q1bc(49:end-1,1) = interfS(1:end-1)*2*dt/Re/dy/dy;
        Q1bc(end,1)      = Q1bc(end,1) + interfS(end)*2*dt/Re/dy/dy;
        Q1bc(end,end)    = Q1bc(end,end) + q1Nlast*2*dt/Re/dy/dy;
        Q1bc             = Q1bc(:);
end