function p = PoissonSolve(p0,rhs,Rp,Ap,Pp,xLen,yLen)

  global dx dy

  %...V-cycles
  rn  = norm(rhs - Ap{1}*p0);
  ic  = 0;
  p   = p0;

  while ( (rn > 1e-8) & (ic < 10) )
    p  = MGVp(1,p,rhs,Rp,Ap,Pp);
    r  = rhs - Ap{1}*p;
    rn = norm(r);
    ic = ic + 1;
    rhs = rhs - sum(r)*(dx/xLen)*(dy/yLen);
    % fprintf('iter %i (p)resnorm %g\n',ic,rn)
  end
