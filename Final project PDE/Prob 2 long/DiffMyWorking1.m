%--------------------------------------
% Driver routine for diffusion
%--------------------------------------

%   clear all
%   close all

  global dx dy
  global dt Re

  %...Reynolds number
  Re = 500;

  %...domain size
  xLen  = 1;
  yLen  = 1;

  N     = 80;
  M     = 64;

  x       = linspace(0,xLen,N+1);
  y       = linspace(0,yLen,M+1);
  xc      = (x(1:end-1)+x(2:end))/2;
  yc      = (y(1:end-1)+y(2:end))/2;
  [yy,xx] = meshgrid(yc,xc);
  dx      = xLen/N;
  dy      = yLen/M;
  dt      = min(dx,dy)/1.5;

  N1 = 64;
  M1 = 32;
  
  N2 = 32;
  M2 = 64;
  
  k = 2;
    
  %...setup the hierarchy of R,A,P
  [Rd1,Pd1]  = getRPd(N1,M1);
  [AAA1,Ad1] = getAd(N1,M1,Rd1,Pd1);
  
  [Rd2,Pd2]  = getRPd(N2,M2);
  [AAA2,Ad2] = getAd(N2,M2,Rd2,Pd2);
    
  %...setup omega1
  q1    = zeros(N1,M1); q1  = q1(:);
 
  q1S =  zeros(N1,1);
  q1N =  ones(N1,1);
  q1W =  zeros(1,M1);
  q1E =  zeros(1,M1);
 
 
  Q1bc       = zeros(N1,M1);
  Q1bc(1,:)  = q1W*k*dt/Re/dx/dx;
  Q1bc(N1,:) = q1E*k*dt/Re/dx/dx;
  Q1bc(:,1)  = Q1bc(:,1)  + q1S*k*dt/Re/dy/dy;
  Q1bc(:,M1) = Q1bc(:,M1) + q1N*k*dt/Re/dy/dy;
  Q1bc       = Q1bc(:);
  
  %...setup omega2
  q2    = zeros(N2,M2); q2  = q2(:);
 
  q2S =  zeros(N2,1);
  q2N =  ones(N2,1);
  q2W =  zeros(1,M2);
  q2E =  zeros(1,M2);
 
 
  Q2bc       = zeros(N2,M2);
  Q2bc(1,:)  = q2W*k*dt/Re/dx/dx;
  Q2bc(N2,:) = q2E*k*dt/Re/dx/dx;
  Q2bc(:,1)  = Q2bc(:,1)  + q2S*k*dt/Re/dy/dy;
  Q2bc(:,M2) = Q2bc(:,M2) + q2N*k*dt/Re/dy/dy;
  Q2bc       = Q2bc(:);
  
     
 sol = zeros(N,M);
 sol(1:48,1:32) = NaN;

 stop = 1e-7;

    interfS = zeros(16,1);
    interfE = zeros(1,32);
    interfW = zeros(1,32);
  for i=1:500
   %...diffusion

    count = 0;
    
    q1old = q1;
    q2old = q2;
    
      while 1
        count = count + 1;
        
        interfSold = interfS;
        interfEold  = interfE;
        interfWold = interfW;
        
        q1    = DiffusionSolve(q1old,Q1bc,AAA1,Rd1,Ad1,Pd1);
        
        sol1 = reshape(q1,N1,M1);
        
        interfW = (sol1(48,:)+sol1(49,:))/2;
        
          q2S =  zeros(N2,1);
          q2N =  ones(N2,1);
          q2W =  zeros(1,M2);
          q2W(33:end) = interfW;
          q2E =  zeros(1,M2);          

          Q2bc       = zeros(N2,M2);
          Q2bc(1,:)  = q2W*k*dt/Re/dx/dx;
          Q2bc(N2,:) = q2E*k*dt/Re/dx/dx;
          Q2bc(:,1)  = Q2bc(:,1)  + q2S*k*dt/Re/dy/dy;
          Q2bc(:,M2) = Q2bc(:,M2) + q2N*k*dt/Re/dy/dy;
          Q2bc       = Q2bc(:);
          
        q2    = DiffusionSolve(q2old,Q2bc,AAA2,Rd2,Ad2,Pd2);
        
        sol2 = reshape(q2,N2,M2);
        
        interfS = (sol2(1:16,32)+sol2(1:16,33))/2;
        interfE = (sol2(16,33:end)+sol2(17,33:end))/2;
       
          q1S =  zeros(N1,1);
          q1S(49:end) = interfS;
          q1N = ones(N1,1);
          q1W = zeros(1,M1);
          q1E = interfE; 

          Q1bc       = zeros(N1,M1);
          Q1bc(1,:)  = q1W*k*dt/Re/dx/dx;
          Q1bc(N1,:) = q1E*k*dt/Re/dx/dx;
          Q1bc(:,1)  = Q1bc(:,1)  + q1S*k*dt/Re/dy/dy;
          Q1bc(:,M1) = Q1bc(:,M1) + q1N*k*dt/Re/dy/dy;
          Q1bc       = Q1bc(:);

        if ( (norm(interfW-interfWold) < stop) && (norm(interfS-interfSold) < stop) && (norm(interfE-interfEold) < stop))
            break;
        end
      end
      display(count);

      
    if (mod(i,500)==0)
      fprintf('time step %i \n',i)
      display(max(abs(interfW-interfWold)));
      display(max(abs(interfS-interfSold)));
      display(max(abs(interfE-interfEold)));
      %...graphics output
%       sol(49:end,:)   = sol2;
%       sol(1:48,33:end) = sol1(1:48,:);
      
      sol(49:end,:)   = sol2;
      sol(1:N1,33:end) = sol1;
      
      figure(1);
      contourf(xx,yy,sol);
      ttl = sprintf('Two domains: Time step %i',i);
%       err = sprintf('maximal error = %.4g',maxerror(i));
%       title([err char(10) ttl])
      colorbar
      title(ttl)
      axis([0 xLen 0 yLen])
      drawnow
%    surf(xx,yy,sol);
%     shading interp;
%     lighting phong;
%     camlight headlight;
%     material shiny;
%     view(0,45)
%     axis([0 xLen 0 yLen 0 1])
%     drawnow
    end
  end
