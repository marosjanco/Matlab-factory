function [qx,qy] = Diff(q,N,M,typ,qS,qN,qW,qE,domain)

global dx dy

  qq = reshape(q,N,M);
  qG = zeros(N+2,M+2);

  qG(2:N+1,2:M+1)   =  qq;
  if (typ=='D')
        qG(1,2:M+1)     = 2*qW-qG(2,2:M+1);
        qG(N+2,2:M+1)   = 2*qE-qG(N+1,2:M+1);
        qG(2:N+1,1)     = 2*qS-qG(2:N+1,2);
        qG(2:N+1,M+2)   = 2*qN-qG(2:N+1,M+1);
        
  elseif (typ=='N')
      
      if (domain==1) 

        qG(1,2:M+1)     = qG(2,2:M+1);
        qG(N+2,2:M+1)   = 2*qE-qG(N+1,2:M+1);
        qG(2:N+1,1)     = qG(2:N+1,2);
        qG(2:N+1,M+2)   = qG(2:N+1,M+1);
      
      elseif (domain==2)
        qG(1,2:33)      = qG(2,2:33);
        qG(1,34:M+1)    = 2*qW(33:end)-qG(2,34:M+1);
        qG(N+2,2:M+1)   = qG(N+1,2:M+1);
        qG(2:N+1,1)     = qG(2:N+1,2);
        qG(2:N+1,M+2)   = qG(2:N+1,M+1);
       else
         fprintf('\nError:\nIncorrect parameter for domain label.\n');
         return; 
      end
  else
      fprintf('\nError:\nIncorrect parameter for Boundary condition label.\n');
      return; 
  end
  qx = (qG(3:N+2,2:M+1)-qG(1:N,2:M+1))/(2*dx);
  qy = (qG(2:N+1,3:M+2)-qG(2:N+1,1:M))/(2*dy);
