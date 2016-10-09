function u = GSd(ilev,u,f,mtimes)
%======================================
% Gauss-Seidel smoother for Diffusion
% equation (matrix-based)
%======================================

  global Rd Ad Pd

  %...Gauss-Seidel split
  AA = Ad{ilev};
  L  = tril(AA);
  U  = triu(AA,1);

  %...iteration
  for k=1:mtimes
    u = L\(f - U*u);
  end
