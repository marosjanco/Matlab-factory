function u = GSp(ilev,u,f,mtimes)
%======================================
% Gauss-Seidel smoother for Poisson
% equation (matrix-based)
%======================================

  global Rp Ap Pp

  %...Gauss-Seidel split
  AA = Ap{ilev};
  L  = tril(AA);
  U  = triu(AA,1);

  %...iteration
  for k=1:mtimes
    u = L\(f - U*u);
  end
