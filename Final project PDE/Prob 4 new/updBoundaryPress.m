function  [P1bc,P2bc] = updBoundaryPress(P1bc,P2bc,p10,N1,M1,p20,N2,M2)

          P1 = reshape(p10,N1,M1);
          P2 = reshape(p20,N2,M2);
          
          P1bc = reshape(P1bc,N1,M1);
          P2bc = reshape(P2bc,N2,M2);
          
          P1bc(49:end,1) = (P2(1:16,32)+P2(1:16,33))/2;
          P1bc(end,:)    = (P2(16,33:end)+P2(17,33:end))/2;
          
          P2bc(1,33:end) = (P1(48,:)+P1(49,:))/2;
          
          P1bc = P1bc(:);
          P2bc = P2bc(:);
end