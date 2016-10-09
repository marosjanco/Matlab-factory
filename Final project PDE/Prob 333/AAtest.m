N1 = 64;
M1 = 32;

N2 = 32;
M2 = 64;


  [Rp1,Pp1] = getRP2p(N1,M1,1);
  Ap1       = getAp(N1,M1,1,Rp1,Pp1);
  ApNiko1   = ApNiko(N1,M1,Rp1,Pp1);
  
  [Rp2,Pp2] = getRP2p(N2,M2,2);
  Ap2       = getAp(N2,M2,2,Rp2,Pp2);
  [ApNiko2,AAA]  = Ap2Niko(N2,M2,Rp2,Pp2);
  
  max(max(abs(full(Ap1{1})-full(ApNiko1{1}))))
  max(max(abs(full(Ap2{1})-full(ApNiko2{1}))))