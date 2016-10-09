%--------------------------------------
% Driver routine for diffusion
%--------------------------------------

  clear all
  close all

  global dx dy
  global dt Re

  %...Reynolds number
  Re = 500;

  %...domain size
  xLen  = 5;
  yLen  = 4;

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
      
  %...setup the hierarchy of R,A,P
  [Rd1,Pd1]  = getRPd(N1,M1);
  [AAA1,Ad1] = getAd(N1,M1,Rd1,Pd1);
  
  [Rd2,Pd2]  = getRPd(N2,M2);
  [AAA2,Ad2] = getAd(N2,M2,Rd2,Pd2);
    
  %...set up omega1
  q1    = zeros(N1,M1); q1  = q1(:);
 
  q1S =  zeros(N1,1);
  q1N =  ones(N1,1);
  q1W =  zeros(1,M1);
  q1E =  zeros(1,M1);
 
 
  Q1bc       = zeros(N1,M1);
  Q1bc(1,:)  = q1W*2*dt/Re/dx/dx;
  Q1bc(N1,:) = q1E*2*dt/Re/dx/dx;
  Q1bc(:,1)  = Q1bc(:,1)  + q1S*2*dt/Re/dy/dy;
  Q1bc(:,M1) = Q1bc(:,M1) + q1N*2*dt/Re/dy/dy;
  Q1bc       = Q1bc(:);
  
  %...set up omega2
  q2    = zeros(N2,M2); q2  = q2(:);
 
  q2S =  zeros(N2,1);
  q2N =  ones(N2,1);
  q2W =  zeros(1,M2);
  q2E =  zeros(1,M2);
 
  Q2bc       = zeros(N2,M2);
  Q2bc(1,:)  = q2W*2*dt/Re/dx/dx;
  Q2bc(N2,:) = q2E*2*dt/Re/dx/dx;
  Q2bc(:,1)  = Q2bc(:,1)  + q2S*2*dt/Re/dy/dy;
  Q2bc(:,M2) = Q2bc(:,M2) + q2N*2*dt/Re/dy/dy;
  Q2bc       = Q2bc(:);
 
  sol = zeros(N,M);
  sol(1:48,1:32) = NaN;

  stop = 1e-5;

  interfSd = zeros(16,1);
  interfEd = zeros(1,32);
  interfWd = zeros(1,32);
  

  
  for i=1:5
   
    count = 0;
    
    q1old = q1;
    q2old = q2;
    
       for j=1:15
    
          count = count + 1;

          sol2    = reshape(q2,N2,M2);
        
          interfS = (sol2(1:16,32)+sol2(1:16,33))/2;
          interfE = (sol2(16,33:end)+sol2(17,33:end))/2;

        
          q1S(49:end) =  interfS;
          q1N =  ones(N1,1);
          q1W =  zeros(1,M1);
          q1E =  interfE;


          Q1bc       = reshape(Q1bc,N1,M1);
          Q1bc(1,:)  = q1W*2*dt/Re/dx/dx;
          Q1bc(N1,:) = q1E*2*dt/Re/dx/dx;
          Q1bc(:,1)  = Q1bc(:,1)  + q1S*2*dt/Re/dy/dy;
          Q1bc(:,M1) = Q1bc(:,M1) + q1N*2*dt/Re/dy/dy;
          Q1bc       = Q1bc(:);
         
          
          q1  = DiffusionSolve(q1old,Q1bc,AAA1,Rd1,Ad1,Pd1);
        

%         Q2bc           = reshape(Q2bc,N2,M2);
%         Q2bc(1,33:end) = interfW*2*dt/Re/dx/dx;
%         Q2bc(1,end)    = Q2bc(1,end) + q2Nfirst*2*dt/Re/dy/dy;
%         Q2bc           = Q2bc(:);
              
          sol1 = reshape(q1,N1,M1);
          interfW = (sol1(48,:)+sol1(49,:))/2;
        
          q2S =  zeros(N2,1);
          q2N =  ones(N2,1);
          q2W(33:end) =  interfW;
          q2E =  zeros(1,M2);
  
          Q2bc       = reshape(Q2bc,N2,M2);
          Q2bc(1,:)  = q2W*2*dt/Re/dx/dx;
          Q2bc(N2,:) = q2E*2*dt/Re/dx/dx;
          Q2bc(:,1)  = Q2bc(:,1)  + q2S*2*dt/Re/dy/dy;
          Q2bc(:,M2) = Q2bc(:,M2) + q2N*2*dt/Re/dy/dy;
          Q2bc       = Q2bc(:);

        q2    = DiffusionSolve(q2old,Q2bc,AAA2,Rd2,Ad2,Pd2);

%         Q1bc             = reshape(Q1bc,N1,M1);
%         Q1bc(end,:)      = interfE*2*dt/Re/dx/dx;
%         Q1bc(49:end-1,1) = interfS(1:end-1)*2*dt/Re/dy/dy;
%         Q1bc(end,1)      = Q1bc(end,1) + interfS(end)*2*dt/Re/dy/dy;
%         Q1bc(end,end)    = Q1bc(end,end) + q1Nlast*2*dt/Re/dy/dy;
%         Q1bc             = Q1bc(:);

          get_error(q1,q2,N1,M1,N2,M2)
       
      end

    
  
    if (mod(i,5)==0)

      fprintf('time step %i \n',i);

            %...putting the sub-solutions together
      sol1 = reshape(q1,N1,M1);
      sol2 = reshape(q2,N2,M2);
      
      sol(49:end,:)   = sol2;
      sol(1:48,33:end) = sol1(1:48,:);    
      
      %...graphics output
      figure(1);
      contourf(xx,yy,sol);
      ttl = sprintf('Two domains: Time step %i',i);
      colorbar
      title(ttl)
      axis([0 xLen 0 yLen])
      drawnow


    end
  end
  