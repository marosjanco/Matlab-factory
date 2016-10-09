function Ap = getAp(N,M,domain,Rp,Pp)
%======================================
% set up a hierarchy of Laplace
% operators (as cell arrays)
%======================================
  global dx dy
  
  %...number of levels and initialization
  kk   = length(Pp);
  Ap   = cell(1,kk+1);

  %...Laplace operator
  AAx      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
  AAx(1,1) = -1/dx/dx;
  AAx(N,N) = -1/dx/dx;

  AAy      = spdiags(ones(M,1)*[1 -2 1]/dy/dy,-1:1,M,M);
  AAy(1,1) = -1/dy/dy;
  AAy(M,M) = -1/dy/dy;
  
  if (domain == 1)
      
          AAx(N,N)    = -3/dx/dx;
          
          AAyDir      = AAy;
          AAyDir(1,1) = -3/dy/dy;
          
          IdNeumx = spdiags([ones(48,1);zeros(16,1)],0,N,N);
          IdDirx  = speye(N)-IdNeumx;

          AAA  = kron(speye(M),AAx) + kron(AAy,IdNeumx) + kron(AAyDir,IdDirx);

  elseif (domain == 2)
          AAxDir      = AAx;
          AAxDir(1,1) = -3/dx/dx;
          
          AAxCorner = AAxDir;
          AAxCorner(N,N) = -3/dx/dx;
          
          AAyCorner = AAy;
          AAyCorner(M,M) = -3/dy/dy;
          
          IdNeumx = spdiags([ones(32,1);zeros(32,1)],0,M,M);
          IdDirx = spdiags([zeros(32,1);ones(31,1);0],0,M,M);
          IdCornerx  = spdiags([zeros(63,1);1],0,M,M);
          
          IdNeumy = spdiags([ones(N-1,1);0],0,N,N);
          IdCornery = speye(N) - IdNeumy;
          
          AAA  = kron(IdNeumx,AAx) + kron(IdDirx,AAxDir)+ kron(IdCornerx,AAxCorner) + kron(AAy,IdNeumy) + kron(AAyCorner,IdCornery);

  else
     fprintf('\nError:\nIncorrect parameter for domain label.\n');
     return; 
  end
      
  %...lower-level operators (Galerkin condition)
  Ap{1} = AAA;
  for i=2:kk+1
    Ap{i} = Rp{i-1}*Ap{i-1}*Pp{i-1};
  end
