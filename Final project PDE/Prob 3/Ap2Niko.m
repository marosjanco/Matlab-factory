function [Ap,AAA] = Ap2Niko(N,M,Rp,Pp)
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
  
  
  AAyND      = spdiags(ones(M,1)*[1 -2 1]/dy/dy,-1:1,M,M);
  AAyND(1,1) = -1/dy/dy;
  AAyND(M,M) = -3/dy/dy;
  
  
  %...Laplace operator Neumann - Dirichlet
  AAxDN      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
  AAxDN(1,1) = -3/dx/dx;
  AAxDN(N,N) = -1/dx/dx;
  
  %...Laplace operator Neumann - Dirichlet
  AAxDD      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
  AAxDD(1,1) = -3/dx/dx;
  AAxDD(N,N) = -3/dx/dx;


  IDx_special = speye(M);IDx_special(1:end-1,1:end-1)=0;
  IDx1 = speye(M);IDx1(33:end,33:end)=0;
  IDx2 = speye(M);IDx2(1:32,1:32)=0;IDx2(end,end)=0;
  IDy = speye(N);IDy(end,end)=0;
  IDy2 = speye(N);IDy2(1:end-1,1:end-1)=0;

  AAA  = kron(IDx_special,AAxDD)+ kron(IDx1,AAxNN) + kron(IDx2,AAxDN) + kron(AAyNN,IDy) + kron(AAyND,IDy2);


  %...lower-level operators (Galerkin condition)
  Ap{1} = AAA;
  for i=2:kk+1
    Ap{i} = Rp{i-1}*Ap{i-1}*Pp{i-1};
  end
