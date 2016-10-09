function qnew = SemiLagrAdvect(u,v,q,qS,qN,qW,qE,N,M)

   global dx dy
   global dt Re

   u = reshape(u,N,M);
   v = reshape(v,N,M);
   q = reshape(q,N,M);

   %...embedding
   qq              = zeros(N+2,M+2);
   qq(2:N+1,2:M+1) = q;

   %...set the ghost values (four edges)
   qq(1,2:M+1)   = 2*qW-qq(2,2:M+1);
   qq(N+2,2:M+1) = 2*qE-qq(N+1,2:M+1);
   qq(2:N+1,1)   = 2*qS-qq(2:N+1,2);
   qq(2:N+1,M+2) = 2*qN-qq(2:N+1,M+1);
   %...set the ghost values (four corners)
   qq(1,1)       = -qq(2,2);
   qq(N+2,1)     = -qq(N+1,2);
   qq(N+2,M+2)   = -qq(N+1,M+1);
   qq(1,M+2)     = -qq(2,M+1);

   q1   = qq(2:N+1,2:M+1);

   q2p  = qq(3:N+2,2:M+1);
   q2m  = qq(1:N,2:M+1);

   q3p  = qq(2:N+1,3:M+2);
   q3m  = qq(2:N+1,1:M);

   q4pp = qq(3:N+2,3:M+2);
   q4mm = qq(1:N,1:M);
   q4pm = qq(3:N+2,1:M);
   q4mp = qq(1:N,3:M+2);

   xi  = -u*dt/dx;
   eta = -v*dt/dy;

   Q2  = q2p.*(xi>0) + q2m.*(xi<0);
   Q3  = q3p.*(eta>0) + q3m.*(eta<0);
   Q4  = q4pp.*((xi>0) & (eta>0)) + q4mm.*((xi<0) & (eta<0)) + ...
         q4pm.*((xi>0) & (eta<0)) + q4mp.*((xi<0) & (eta>0));

   qnew = (1-abs(xi)).*(1-abs(eta)).*q1 + ...
          abs(xi).*(1-abs(eta)).*Q2 + ...
	  abs(eta).*(1-abs(xi)).*Q3 + ...
	  abs(xi).*abs(eta).*Q4;


   qnew = qnew(:);
