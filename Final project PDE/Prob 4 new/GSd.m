function u = GSd(AA,u,f,mtimes)
%======================================
% Gauss-Seidel smoother for Diffusion
% equation (matrix-based)
%======================================

  %...Gauss-Seidel split
  L  = tril(AA);
  U  = triu(AA,1);

  %...iteration
  for k=1:mtimes
    u = L\(f - U*u);
  end