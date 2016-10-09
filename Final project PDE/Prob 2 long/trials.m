  N = 4;
  M = 3;
  
  vS = (1:N)';%zeros(N2,1);
  vN = (N:-1:1)';%zeros(N2,1);
  vW = 1:M;%zeros(1,M2);
  vE = M:-1:1;%zeros(1,M2);

  Vbc      = zeros(N,M);
  Vbc(1,:) = vW;
  Vbc(N,:) = vE;
  Vbc(:,1) = vS;
  Vbc(:,M) = vN;
  Vbc
  Vbc = Vbc(:);
  Vbc