%--------------------------------------
% Driver routine for diffusion
%--------------------------------------
 
  clear all
  close all
 
  global dx dy
  global dt Re
 
  %...Reynolds number
  Re = 400;
 
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
 
  sol1 = reshape(q1,N1,M1);
  sol2 = reshape(q2,N2,M2);
  
  interfS = zeros(16,1);
  interfE = zeros(1,32);
  interfW = zeros(1,32);
  
  stop = 1e-6;       %stopping criterion for multiplicative Schwarz
  TimeEnd = 2000;    %Given in original DriverMain.m
  
  % For Residual Plots
  plotStep = 50;
  Schwarz_maxiter = 4; %maximal number for multiplicative Schwarz to be considered for plotting residuals
  Timesteps = plotStep:plotStep:TimeEnd;
  Residuals = zeros(length(Timesteps),Schwarz_maxiter);
  Residuals200 = zeros(200,Schwarz_maxiter);
  
  for i=1:TimeEnd
    
    Residual  = zeros(1,Schwarz_maxiter); % For Residual Plots
    count = 0;
    
    q1old = sol1(:);
    q2old = sol2(:);
    
      while 1
        count = count + 1;
        
        interfSold = interfS;
        interfEold = interfE;
        interfWold = interfW;
    
        q1    = DiffusionSolve(q1old,Q1bc,AAA1,Rd1,Ad1,Pd1);
        [Q2bc, interfW] = setBoundary2d(q1,N1,M1,Q2bc,N2,M2,q2N(1));
        q2    = DiffusionSolve(q2old,Q2bc,AAA2,Rd2,Ad2,Pd2);
        [Q1bc, interfS, interfE] = setBoundary1d(q2,N2,M2,Q1bc,N1,M1,q1N(end));
       
        Res = norm([(interfW-interfWold),(interfS-interfSold)',(interfE-interfEold)]);
            
        % For Residual Plots
        if (count<=Schwarz_maxiter)
        Residual(count) = Res;
        end
        
        if (Res < stop)
            break;
        end
      end
      
      %unify OL
      [q1,q2] = unify(q1,N1,M1,q2,N2,M2);

      display(count);
      
      sol2 = reshape(q2,N2,M2);
      sol1 = reshape(q1,N1,M1);
      
      sol1(57:end,:) = sol2(9:16,33:end);
      sol2(1:8,33:end) = sol1(49:56,:);
      
    % For Residual Plots
    if i<=200
        Residuals200(i,:) = Residual;
    end
    
    if (mod(i,plotStep)==0)
      fprintf('time step %i \n',i);
      
      % For Residual Plots
      Residuals(fix(i/plotStep),:) = Residual;
      
      %...putting the sub-solutions together    
      sol(49:end,:)   = sol2;
      sol(1:48,33:end) = sol1(1:48,:);    
      
      %...graphics output
      figure(1);
      contourf(xx,yy,sol);
      colorbar
      title(sprintf('Composite solution at time-step %i',i))
      axis([0 xLen 0 yLen])
      drawnow
    end
  end
  
  % For Residual Plots
  
  % Residuals evolving over time-steps 50*k, where k = 1,2,3,...
  [~, yloc] = ind2sub(size(max(Residuals)),find(max(Residuals)==0));
  NumRes = yloc(1)-1;
  Residuals(Residuals==0) = NaN;
  
  figure(2)
  hold on
  for i = 1:NumRes
      subplot(NumRes,1,i)
      plot(Timesteps,Residuals(:,i))
      if i == 1
          title(sprintf('   Residuals evolving over time-steps %d*k, k = 1,2,3,...',plotStep))
      end
      ylabel('Residual')
      legend(sprintf('Residuals after %d. exchange in overlapping region',i),'bestoutside')
  end
  xlabel('Time-step')
  hold off
 
 % Residuals evolving over the first 200 time-steps
  [~, yloc] = ind2sub(size(max(Residuals200)),find(max(Residuals200)==0));
  NumRes = yloc(1)-1;
  Residuals200(Residuals200==0) = NaN;
  
  figure(3)
  hold on
  for i = 1:NumRes
      subplot(NumRes,1,i)
      plot(1:200,Residuals200(:,i))
      if i == 1
            title('Residuals evolving over the first 200 time-steps')
      end
      ylabel('Residual')
      lgnd = sprintf('Residuals after %d. exchange in overlapping region',i);
      legend(lgnd,'bestoutside')
  end
 xlabel('Time-step')
 hold off

