function [u1new,u2new] = unify(u1new,N1,M1,u2new,N2,M2)

            soltemp2 = reshape(u2new,N2,M2);
            soltemp1 = reshape(u1new,N1,M1);
          
            soltemp1(57:end,:) = soltemp2(9:16,33:end);
            soltemp2(1:8,33:end) = soltemp1(49:56,:);
          
            u1new = soltemp1(:);
            u2new = soltemp2(:);
end