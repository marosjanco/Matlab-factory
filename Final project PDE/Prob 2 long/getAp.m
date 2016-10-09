function getAp(N,M)
%======================================
% set up a hierarchy of Laplace
% operators (as cell arrays)
%======================================

  global Rp Ap Pp
  global xLen yLen

  %...number of levels and initialization
  kk   = length(Pp);
  dx   = xLen/N;
  dy   = yLen/M;
  Ap   = cell(1,kk+1);

  %...Laplace operator
  AAx      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
  AAx(1,1) = -1/dx/dx;
  AAx(N,N) = -1/dx/dx;

  AAy      = spdiags(ones(M,1)*[1 -2 1]/dy/dy,-1:1,M,M);
  AAy(1,1) = -1/dy/dy;
  AAy(M,M) = -1/dy/dy;

  AAA  = kron(speye(M),AAx) + kron(AAy,speye(N));

  %...lower-level operators (Galerkin condition)
  Ap{1} = AAA;
  for i=2:kk+1
    Ap{i} = Rp{i-1}*Ap{i-1}*Pp{i-1};
  end
