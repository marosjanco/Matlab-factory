function u = MGVp(ilev,u,rhs,Rp,Ap,Pp)
%======================================
% multigrid V-cycle (recursive)
%======================================

  ilevmin = length(Rp);
  mu = 5;

  %...lowest level
  if (ilev==ilevmin+1)
    u  = GSp(Ap{ilev},u,rhs,mu);
  else
    %...presmoothing
    u  = GSp(Ap{ilev},u,rhs,mu);

    %...residual and restriction
    r  = rhs - Ap{ilev}*u;
    rr = Rp{ilev}*r;

    %...call MGV
    ee = zeros(size(rr));
    ee = MGVp(ilev+1,ee,rr,Rp,Ap,Pp);

    %...prolong and add
    u  = u + Pp{ilev}*ee;

    %...postsmoothing
    u  = GSp(Ap{ilev},u,rhs,mu);
  end
