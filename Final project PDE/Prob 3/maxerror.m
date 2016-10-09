function err = maxerror(q1,N1,M1,q2,N2,M2)

    q1_2D = reshape(q1,N1,M1);
    q2_2D = reshape(q2,N2,M2); 
    
    err = max(max(abs(q1_2D(49:end,1:32)-q2_2D(1:16,33:end))));

end