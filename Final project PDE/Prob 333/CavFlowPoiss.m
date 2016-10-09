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
      
  x1Len = xLen*N1/N;
  y1Len = yLen*M1/M;
  x2Len = xLen*N2/N;
  y2Len = yLen*M2/M;
  
  %...set up the hierarchy of R,A,P
  
  [Rp1,Pp1] = getRPp(N1,M1,1);
  [Rp2,Pp2] = getRPp(N2,M2,2);
  
  Ap1       = getAp(N1,M1,1,Rp1,Pp1);
  Ap2       = getAp(N2,M2,2,Rp2,Pp2);
  
  %...set up omega1
  q1    = zeros(N1,M1); q1  = q1(:);
  Q1bc = zeros(N1,M1); Q1bc = Q1bc(:);
  
  %...set up omega2
  q2    = zeros(N2,M2); q2  = q2(:);
  Q2bc = zeros(N2,M2); Q2bc = Q2bc(:);
  
  f1 = ones(N1,M1);
  f1 = f1(:);
  
  f2 = ones(N2,M2);
  f2 = f2(:);
  
  
  sol = zeros(N,M);
  sol(1:48,1:32) = NaN;

  stopp = 1e-6;

  interfS = zeros(16,1);
  interfE = zeros(1,32);
  interfW = zeros(1,32);
  
  % For plotting residuals
  Residuals = zeros(20,1);
  
count = 0;

q1old = q1;
q2old = q2;

err = NaN*zeros(20,1);

while 1
    count = count + 1;

    interfSold = interfS;
    interfEold = interfE;
    interfWold = interfW;
    
    q1    = PoissonSolve(q1old,f1-Q1bc,Rp1,Ap1,Pp1,x1Len,y1Len);
    [Q2bc, interfW] = setBoundary2p(q1,N1,M1,Q2bc,N2,M2);
    q2    = PoissonSolve(q2old,f2-Q2bc,Rp2,Ap2,Pp2,x2Len,y2Len);
    [Q1bc, interfS, interfE] = setBoundary1p(q2,N2,M2,Q1bc,N1,M1);

    err(count) = maxerror(q1,N1,M1,q2,N2,M2);
    
    Res = norm([(interfW-interfWold),(interfS-interfSold)',(interfE-interfEold)]);

    % For plotting residuals
    Residuals(count) = Res;

    if Res < stopp
        break;
    end
end

%...putting the sub-solutions together

  sol1 = reshape(q1,N1,M1);
  sol2 = reshape(q2,N2,M2);

  %...graphics output
 
  figure(1)
  
  subplot(2,2,1)
  sol(49:end,:)   = sol2;
  sol(1:48,33:end) = sol1(1:48,:);  
  soll2 = sol;
  contourf(xx,yy,sol);
  colorbar
  title(['Contour-plot soution to' char(10) 'Poisson Equation with Omega 2 on top'])
  axis([0 xLen 0 yLen])
  axis equal

  subplot(2,2,3)
  surf(xx,yy,sol);
  title(['Surface-plot soution to' char(10) 'Poisson Equation with Omega 2 on top'])
  axis tight
  view(-35,35)
    
  subplot(2,2,2)
  sol(49:end,:)   = sol2;
  sol(1:N1,33:end) = sol1;
  
  contourf(xx,yy,sol);
  colorbar
  title(['Contour-plot soution to' char(10) 'Poisson Equation with Omega 1 on top'])
  axis([0 xLen 0 yLen])
  axis equal
  
  subplot(2,2,4)
  surf(xx,yy,sol);
  title(['Surface-plot soution to' char(10) 'Poisson Equation with Omega 1 on top'])
  axis tight
  view(-35,35)

 % For plotting residuals
 figure(2)
 subplot(1,2,1)
 plot(1:count,Residuals(1:count),'-xg');
 title(['Maximal residual over the interfaces' char(10) '(i.e. the west, east, and south boundaries of the Overlap region)']);

 ylabel('Residual')
 xlabel('Iteration number in the while loop of Schwarz algorithm')
 
 subplot(1,2,2)
 maxerr = max(max(abs(soll2-sol)));
 plot(1:length(err),err,'-o');
 title(['Maximal absolute difference (error) in the Overlap region' char(10) sprintf('The maximal error over all iteration is %g',maxerr)]);
 ylabel('Error')
 xlabel('Iteration number in the while loop of Schwarz algorithm')
 
 display(err)
 display(Residuals)
 display(count)
 