function u = MGVp(ilev,u,rhs)
%======================================
% multigrid V-cycle (recursive)
%======================================

  global Rp Ap Pp
  global ilevmin

  mu = 5;

  %...lowest level
  if (ilev==ilevmin+1)
    u  = GSp(ilev,u,rhs,mu);
  else
    %...presmoothing
    u  = GSp(ilev,u,rhs,mu);

    %...residual and restriction
    r  = rhs - Ap{ilev}*u;
    rr = Rp{ilev}*r;

    %...call MGV
    ee = zeros(size(rr));
    ee = MGVp(ilev+1,ee,rr);

    %...prolong and add
    u  = u + Pp{ilev}*ee;

    %...postsmoothing
    u  = GSp(ilev,u,rhs,mu);
  end
