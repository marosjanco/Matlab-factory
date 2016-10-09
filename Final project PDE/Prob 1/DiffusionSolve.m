function unew = DiffusionSolve(u,Ubc)

  global Rd Ad Pd
  global dt Re

  %...setup right-hand side (for Crank-Nicolson)
  rhs = u + dt/2/Re*Ad{1}*u + Ubc;

  rn  = norm(rhs - Ad{1}*u);
  ic  = 0;
  while ( (rn > 1e-10) & (ic < 10) )
    u  = MGVd(1,u,rhs);
    rn = norm(rhs - Ad{1}*u);
    ic = ic + 1;
  end
  %fprintf('iter %i (d)resnorm %g\n',ic,rn)

  unew = u;
