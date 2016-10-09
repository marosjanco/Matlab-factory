     N = 8;
     M = 4;
 
  %...number of levels and initialization
  kk      = log(N)/log(2)-1;
  ll      = log(M)/log(2)-1;
  kl      = min(kk,ll);

  %...set up prolongation
  NN = N;
  MM = M;
  
    PPx = sparse(2*NN,NN);
    for j=1:NN-1
      PPx(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
    end
    PPx(1,1)     = 1;
    PPx(2*NN,NN) = 0.5;
    
    PPy = sparse(2*MM,MM);
    for j=1:MM-1
      PPy(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
    end
    PPy(1,1)     = 0.5;
    PPy(2*MM,MM) = 1;

    xh = full(PPx);
    yh =  full(PPy);
    
    PPP = kron(PPy,PPx);
    PPh = full(PPP);
  