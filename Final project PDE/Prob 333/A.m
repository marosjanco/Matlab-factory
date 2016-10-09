  N = 32;
  M = 64;

  %...number of levels and initialization
  dx   = 1;
  dy   = 1;

  %...Laplace operator
  AAx      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
  AAx(1,1) = -1/dx/dx;
  AAx(N,N) = -1/dx/dx;
  
  AAy      = spdiags(ones(M,1)*[1 -2 1]/dy/dy,-1:1,M,M);
  AAy(1,1) = -1/dy/dy;
  AAy(M,M) = -1/dy/dy;
% 
%   Id1 = zeros(M,M);
%   Id1(1,1) = 1;
%   full(AAx);
%   full(kron(Id1,AAx));
%   
%   Id2 = zeros(N,N);
%   Id2(1,1) = 1;
%   full(AAy);
%   full(kron(AAy,Id2));
%   

  AAxDirichlet = AAx;
  AAxDirichlet(1,1)     = -3/dx/dx;

  IdNeumann    = sparse(M,M);
  for i=1:32
      IdNeumann(i,i) = 1;
  end

  IdDirichlet  = sparse(M,M);
  for i=33:M
      IdDirichlet(i,i) = 1;
  end

AAA  = kron(IdNeumann,AAx) + kron(IdDirichlet,AAxDirichlet) + kron(AAy,speye(N));

%           
%   AAA  = kron(speye(M),AAx) + kron(AAy,speye(N));
%   
%   AAxf = full(AAx)
%   AAyf = full(AAy)
%   
%   AAAf = full(AAA)