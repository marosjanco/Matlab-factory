%--------------------------------------
% Driver routine for Poisson
%--------------------------------------

  clear all
  close all

  global dx dy
  global dt
  global xLen yLen
  
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
      
  %...set up the hierarchy of R,A,P
  [Rp1,Pp1] = getRPp(N1,M1,1);
  Ap1       = getAp(N1,M1,1,Rp1,Pp1);
  
  [Rp2,Pp2] = getRPp(N2,M2,2);
  Ap2       = getAp(N2,M2,2,Rp2,Pp2);
  
  %...set up omega1
  q1    = zeros(N1,M1); q1  = q1(:);
 
  q1S =  zeros(N1,1);
  q1N =  ones(N1,1);
  q1W =  ones(1,M1);
  q1E =  ones(1,M1);
  
  Q1bc       = zeros(N1,M1);
%   Q1bc(1,:)  = q1W/dx;
%   Q1bc(N1,:) = q1E*2/dx/dx;
%   
%   Q1bc(1:48,1)  = Q1bc(1:48,1)  + q1S(1:48)/dy;
%   Q1bc(49:N1,1) = Q1bc(49:N1,1)  + q1S(49:N1)*2/dy/dy;
%   Q1bc(:,M1) = Q1bc(:,M1) + q1N/dy;
  Q1bc       = Q1bc(:);
  
  %...set up omega2
  q2    = zeros(N2,M2); q2  = q2(:);
  q2S =  zeros(N2,1);
  q2N =  ones(N2,1);
  q2W =  ones(1,M2);
  q2E =  ones(1,M2);
 
  Q2bc       = ones(N2,M2);
%   Q2bc(1,1:32)  = q2W(1:32)/dx;
%   Q2bc(1,33:N2)  = q2W(33:N2)*2/dx/dx;
%   Q2bc(N2,:) = q2E/dx;
%   Q2bc(:,1)  = Q2bc(:,1)  + q2S/dy;
%   Q2bc(:,M2) = Q2bc(:,M2) + q2N/dy;
 
  Q2bc       = Q2bc(:);
  sol = zeros(N,M);
  sol(1:48,1:32) = NaN;

  stopp = 1e-6;

  interfSp = zeros(16,1);
  interfEp = zeros(1,32);
  interfWp = zeros(1,32);
  
  % For plotting residuals
  ResTimesteps = [1,5,10,30,50];
  MaxIter = 22;
  NumRes = length(ResTimesteps);
  Residuals = zeros(NumRes,MaxIter);
  iRes = 1;
  
      new1 = Q1bc;
      new2 = Q2bc;
      
  for i=1:50

    count = 0;
    
    q1old = q1;
    q2old = q2;
    
      while 1
        count = count + 1;
        
        interfSpold = interfSp;
        interfEpold = interfEp;
        interfWpold = interfWp;
        
        
        q1    = PoissonSolve(q1old,-new1,Rp1,Ap1,Pp1);
        [new2, interfWp] = setBoundary2p(q1,N1,M1,Q2bc,N2,M2);
        q2    = PoissonSolve(q2old,-new2,Rp2,Ap2,Pp2);
        [new1, interfSp, interfEp] = setBoundary1p(q2,N2,M2,Q1bc,N1,M1);

        Res = max([norm(interfWp-interfWpold), norm(interfSp-interfSpold), norm(interfEp-interfEpold)]);

        % For plotting residuals
        if i == ResTimesteps(iRes)
            Residuals(iRes,count) = Res;
        end
        
        if Res < stopp
            break;
        end
      end
      
          pp1 = reshape(q1,N1,M1);
          pp2 = reshape(q2,N2,M2);
          mx1 = pp1(49:end,:);
          mx2 = pp2(1:16,33:end);
          subplot(2,2,1);
          contourf(mx1);
          subplot(2,2,2);
          contourf(mx2);pause(3);
          
        % For plotting residuals
        if (i == ResTimesteps(iRes) && iRes < NumRes)
            iRes = iRes + 1;
        end
      
    if max(i==ResTimesteps)
      fprintf('time step %i \n',i)
      
      
      %...putting the sub-solutions together
      
      sol1 = reshape(q1,N1,M1);
      sol2 = reshape(q2,N2,M2);
      
      subplot(2,2,3)
      sol(49:end,:)   = sol2;
      sol(1:48,33:end) = sol1(1:48,:);  
      contourf(xx,yy,sol);
      colorbar
      title(sprintf('Omega 2 on top; Time-step %i',i))
      axis([0 xLen 0 yLen])
      axis equal
      
      subplot(2,2,4)
      sol(49:end,:)   = sol2;
      sol(1:N1,33:end) = sol1;
           
      %...graphics output
      contourf(xx,yy,sol); pause(4);
      title(sprintf('Omega 1 on top; Time-step %i',i))
      colorbar
      axis([0 xLen 0 yLen])
      axis equal
      drawnow
    end
  end

    % For plotting residuals
  
  Residuals(Residuals==0) = NaN;
  
  figure(2)
  hold on
  Legend = cell(1,NumRes);
  
  for i = 1:NumRes
      plot(1:MaxIter,Residuals(i,:))
      Legend{i} = horzcat('Residuals at time-step ',num2str(ResTimesteps(i)));
  end
 
 title('Residuals evolving at fixed time-steps = 1,5,10,30,50')
 ylabel('Residual')
 xlabel('Iteration number in the while loop of Schwarz algorithm')
 legend(Legend,'bestoutside')
 hold off
