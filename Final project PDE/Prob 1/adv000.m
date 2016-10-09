%--------------------------------------
% Driver routine for ADVECTION (with error-checking)
%--------------------------------------

  clear all
  close all
  
  global dx dy
  global N M
  global dt

  %...Setting advection velocities (now constants for test purposes)
  c1 = 1;    % in x direction
  c2 = -0.05; % in y direction
  
  %...Domain omega_1
  
  N1 = 48;
  M1 = 32;
  
  ax1 = c1*ones(N1*M1,1);
  ay1 = c2*ones(N1*M1,1);

  q1S =  zeros(N1,1);
  q1N =  zeros(N1,1);
  q1W =  zeros(1,M1);
  q1E =  zeros(1,M1);%  
  
  %...Domain omega_1

  N2 = 32;
  M2 = 64;
  
  ax2 = c1*ones(N2*M2,1);
  ay2 = c2*ones(N2*M2,1);

  q2S =  zeros(N2,1);
  q2N =  zeros(N2,1);
  q2W =  zeros(1,M2);
  q2E =  zeros(1,M2);

  %...domain size
  xLen  = 1;
  yLen  = 0.8;
  
  N = N1 + N2;
  M = M2;
  
  x       = linspace(0,xLen,N+1);
  y       = linspace(0,yLen,M+1);
  xc      = (x(1:end-1)+x(2:end))/2;
  yc      = (y(1:end-1)+y(2:end))/2;
  [yy,xx] = meshgrid(yc,xc);
  
  dx      = xLen/N;
  dy      = yLen/M;
  dt      = min(dx,dy)/1.5;
  
  
  f_init = @(x,y) exp(-100.*((x-0.3).^2+(y-0.6).^2));
  sol = f_init(xx,yy);
  sol(1:N1,1:M2-M1) = NaN;
  
  sol1 = sol(1:N1,M2-M1+1:end);
  sol2 = sol(N1+1:end,:);

  q1 = sol1(:);
  q2 = sol2(:);
  
  % Setting up true advection for omega1 union omega2 (without domain decomposion) 
  soltrue = sol(:,M2-M1+1:end);
  qtrue   = soltrue(:);
  soltrue = zeros(size(sol));
  soltrue(:,M2-M1+1:end) = sol(:,M2-M1+1:end);
  soltrue(:,1:M2-M1) = NaN;
  
  c1qtrue = c1*ones(size(qtrue));
  c2qtrue = c2*ones(size(qtrue));

  qtrueS =  zeros(N,1);
  qtrueN =  zeros(N,1);
  qtrueW =  zeros(1,M1);
  qtrueE =  zeros(1,M1);
  
  maxerror = zeros(90,1);
  for i=1:60

    %...Updating the boundaries between omega_1 and omega_2
    q1E  = (sol1(N1,:)+sol2(1,M2-M1+1:end))/2;
    q2W(M2-M1+1:end) = q1E;
    
    %...Advecting omega_1 and omega_2 separately
    q1   = SemiLagrAdvect(ax1,ay1,q1,q1S,q1N,q1W,q1E,N1,M1);   
    q2    = SemiLagrAdvect(ax2,ay2,q2,q2S,q2N,q2W,q2E,N2,M2);
    
    %...Extracting sub-solutions
    sol1 = reshape(q1,N1,M1);
    sol2  = reshape(q2,N2,M2);
      
    %...Putting the sub-solutions together
    sol(1:N1,M2-M1+1:end) = sol1;
    sol(N1+1:end,:) = sol2;
    
    qtrue    = SemiLagrAdvect(c1qtrue,c2qtrue,qtrue,qtrueS,qtrueN,qtrueW,qtrueE,N,M1);
    soltrue(:,M2-M1+1:end) = reshape(qtrue,N,M1);

    error = abs(sol - soltrue);
    maxerror(i) = max(max(error));

    if (i==1 || (mod(i,10)==0))
      fprintf('time step %i \n',i)
      
      figure(1)
      subplot(2,1,1)   
      contourf(xx,yy,soltrue);
      ttl = sprintf('One domain: Time step %i',i);
      title(ttl)
      colorbar
      xlabel('x')
      ylabel('y')
      axis image;
      
      subplot(2,1,2)
      contourf(xx,yy,sol);
      ttl = sprintf('Two domains: Time step %i',i);
      err = sprintf('maximal error = %.4g',maxerror(i));
      title([err char(10) ttl])
      colorbar
      xlabel('x')
      ylabel('y')
      axis image;
      drawnow
    end
  end
  display(max(maxerror),'maximal errror');