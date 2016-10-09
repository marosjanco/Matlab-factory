function u = MGVd(ilev,u,rhs)
%======================================
% multigrid V-cycle (recursive)
%======================================

  global Rd Ad Pd
  global ilevmin

  mu = 5;

  %...lowest level
  if (ilev==ilevmin+1)
    u  = Ad{ilev}\rhs;
  else
    %...presmoothing
    u  = GSd(ilev,u,rhs,mu);

    %...residual and restriction
    r  = rhs - Ad{ilev}*u;
    rr = Rd{ilev}*r;

    %...call MGV
    ee = zeros(size(rr));
    ee = MGVd(ilev+1,ee,rr);

    %...prolong and add
    u  = u + Pd{ilev}*ee;

    %...postsmoothing
    u  = GSd(ilev,u,rhs,mu);
  end
