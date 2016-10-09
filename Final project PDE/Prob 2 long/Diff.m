function [qx,qy] = Diff(q,N,M,typ,qS,qN,qW,qE)

  global xLen yLen

  dx = xLen/N;
  dy = yLen/M;

  qq = reshape(q,N,M);
  qG = zeros(N+2,M+2);

  qG(2:N+1,2:M+1)   =  qq;
  if (typ=='D')
    qG(1,2:M+1)     = 2*qE-qG(2,2:M+1);
    qG(N+2,2:M+1)   = 2*qW-qG(N+1,2:M+1);
    qG(2:N+1,1)     = 2*qS-qG(2:N+1,2);
    qG(2:N+1,M+2)   = 2*qN-qG(2:N+1,M+1);
  elseif (typ=='N')
    qG(1,2:M+1)     = qG(2,2:M+1);
    qG(N+2,2:M+1)   = qG(N+1,2:M+1);
    qG(2:N+1,1)     = qG(2:N+1,2);
    qG(2:N+1,M+2)   = qG(2:N+1,M+1);
  end

  qx = (qG(3:N+2,2:M+1)-qG(1:N,2:M+1))/(2*dx);
  qy = (qG(2:N+1,3:M+2)-qG(2:N+1,1:M))/(2*dy);
