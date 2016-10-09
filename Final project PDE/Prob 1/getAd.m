function getAd(N,M)
%======================================
% set up a hierarchy of Laplace
% operators (as cell arrays)
%======================================

  global Rd Ad Pd
  global xLen yLen
  global dt Re

  %...number of levels and initialization
  kk   = length(Pd);
  dx   = xLen/N;
  dy   = yLen/M;
  A    = cell(1,kk+1);

  %...Laplace operator
  AAx      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
  AAx(1,1) = -3/dx/dx;
  AAx(N,N) = -3/dx/dx;

  AAy      = spdiags(ones(M,1)*[1 -2 1]/dy/dy,-1:1,M,M);
  AAy(1,1) = -3/dy/dy;
  AAy(M,M) = -3/dy/dy;

  AAA  = kron(speye(M),AAx) + kron(AAy,speye(N));
  DDD  = speye(N*M) - dt/2/Re*AAA;

  %...lower-level operators (Galerkin condition)
  Ad{1} = DDD;
  for i=2:kk+1
    Ad{i} = Rd{i-1}*Ad{i-1}*Pd{i-1};
  end
