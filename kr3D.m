function A = kr3D(n)
    NN = n*n;
    NNN = NN*n;
    
    A = zeros(NNN,NNN);
    
    for i = 1 : NNN
        A(i,i) = -8;
    end
    
    for i = 1 : n
        for j = 1 : NN-n
            A(NN*(i-1)+j,NN*(i-1)+j+n) = 2;
            A(NN*(i-1)+j+n,NN*(i-1)+j) = 2;
        end
    end
    
    for i = 1 : NN
        for j = 1:n-1
            A((i-1)*n +j+1,(i-1)*n +j) = 1;
            A((i-1)*n +j,(i-1)*n +j+1) = 1;
        end
    end
    
    
    for j = NN + 1 : NNN;
        A(j-NN,j) = 1;
        A(j,j-NN) = 1;
    end
        
end











































