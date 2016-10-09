function C = krc3D(n)
    NN = n*n;
    NNN = NN*n;
    
    C = zeros(NNN,2*NN+1);
    
    for i = 1 : NNN
        C(i,NN+1) = -8;
    end
    
    for i = 1 : n
        for j = 1 : NN-n
            C(NN*(i-1)+j+n, NN + 1 - n) = 2;
            C(NN*(i-1)+j,NN + 1 + n) = 2;
        end
    end
    
    for i = 1 : NN
        for j = 1:n-1
            C((i-1)*n +j+1,NN) = 1;
            C((i-1)*n +j,NN+2) = 1;
        end
    end
    
    
    for j = NN + 1 : NNN;
        C(j-NN,2*NN+1) = 1;
        C(j,1) = 1;
    end
        
end











































