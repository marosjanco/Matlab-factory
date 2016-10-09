function y = kr(n1,B)
    
    N = n1*B;
    y = zeros(N,N);
    
    for i = 1 : N
        y(i,i) = -4;
    end
    
    for i = 1 : B
        for j = 1:n1-1
            y((i-1)*n1 + j+1,(i-1)*n1 +j) = 1;
            y((i-1)*n1 +j,(i-1)*n1 +j+1) = 1;
        end
    end
    for j = B + 1 : N;
        y(j-B,j) = 1;
        y(j,j-B) = 1;
    end
        
end