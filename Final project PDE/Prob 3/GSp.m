function u = GSp(AA,u,f,mtimes)
%======================================
% Gauss-Seidel smoother for Poisson
% equation (matrix-based)
%======================================

  %...Gauss-Seidel split
  L  = tril(AA);
  U  = triu(AA,1);

  %...iteration
  for k=1:mtimes
    u = L\(f - U*u);
  end