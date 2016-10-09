function [Ap,AAA] = ApNiko(N,M,Rp,Pp)
%======================================
% set up a hierarchy of Laplace
% operators (as cell arrays)
%======================================
  global dx dy
  %...number of levels and initialization
  kk   = length(Pp);
  Ap   = cell(1,kk+1);
  %...Laplace operator Neumann-Neumman
  AAxNN      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
  AAxNN(1,1) = -1/dx/dx;
  AAxNN(N,N) = -1/dx/dx;
  AAyNN      = spdiags(ones(M,1)*[1 -2 1]/dy/dy,-1:1,M,M);
  AAyNN(1,1) = -1/dy/dy;
  AAyNN(M,M) = -1/dy/dy;
  
  %...Laplace operator Neu - Dir
  AAxND      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
  AAxND(1,1) = -1/dx/dx;
  AAxND(N,N) = -3/dx/dx;
  
  %...Laplace operator Dir - Neu
  AAyDN      = spdiags(ones(M,1)*[1 -2 1]/dy/dy,-1:1,M,M);
  AAyDN(1,1) = -3/dy/dy;
  AAyDN(M,M) = -1/dy/dy;
 
 % IDx1 = speye(M);
  IDx2 = speye(M);
  IDy1 = speye(N);IDy1(49:end,49:end)=0;
  IDy2 = speye(N);IDy2(1:48,1:48)=0;
  
  AAA  =  kron(IDx2,AAxND) + kron(AAyNN,IDy1) + kron(AAyDN,IDy2);%kron(IDx1,AAxNN) +
  %...lower-level operators (Galerkin condition)
  Ap{1} = AAA;
  for i=2:kk+1
    Ap{i} = Rp{i-1}*Ap{i-1}*Pp{i-1};
  end
