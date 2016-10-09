function y = krb(n1,B) % n = n1 * n2 = n1 * B
    
    N = n1*B;
    y = zeros(N,2*B+1);
    
    for i = 1 : N
        y(i,B+1) = -4;
    end
    
    for i = 1 : B
        for j = 1:n1-1
            y((i-1)*n1 + j+1,B) = 1;
            y((i-1)*n1 +j,B+2) = 1;
        end
    end
    
    for j = B + 1 : N;
        y(j,1) = 1;
        y(j-B,2*B+1) = 1;
    end
        
end