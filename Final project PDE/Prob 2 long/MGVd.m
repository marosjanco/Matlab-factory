function u = MGVd(ilev,u,rhs,Rd,Ad,Pd)
%======================================
% multigrid V-cycle (recursive)
%======================================
  
  ilevmin = length(Rd);
  mu = 5;

  %...lowest level
  if (ilev==ilevmin+1)
    u  = Ad{ilev}\rhs;
  else
    %...presmoothing
    u  = GSd(Ad{ilev},u,rhs,mu);

    %...residual and restriction
    r  = rhs - Ad{ilev}*u;
    rr = Rd{ilev}*r;

    %...call MGV
    ee = zeros(size(rr));
    ee = MGVd(ilev+1,ee,rr,Rd,Ad,Pd);

    %...prolong and add
    u  = u + Pd{ilev}*ee;

    %...postsmoothing
    u  = GSd(Ad{ilev},u,rhs,mu);
  end
