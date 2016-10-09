%--------------------------------------
% Driver routine for Poisson
%--------------------------------------

  clear all
  close all

  global dx dy
  global dt
  global xLen yLen
  
  %...domain size
  xLen  = 80/64;
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
  [Rp1,Pp1] = getRP2p(N1,M1,1);
  Ap1       = getAp(N1,M1,1,Rp1,Pp1);
  
  [Rp2,Pp2] = getRP2p(N2,M2,2);
  Ap2       = getAp(N2,M2,2,Rp2,Pp2);
  
  %...set up omega1
  q1    = zeros(N1,M1); q1  = q1(:);
 
  q1S =  zeros(N1,1);
  q1N =  zeros(N1,1);
  q1W =  zeros(1,M1);
  q1E =  zeros(1,M1);
  
  Q1bc       = zeros(N1,M1);
  Q1bc(1,:)  = q1W/dx;
  Q1bc(N1,:) = q1E*2/dx/dx;
  Q1bc(1:48,1)  = Q1bc(1:48,1)  + q1S(1:48)/dy;
  Q1bc(49:N1,1) = Q1bc(49:N1,1)  + q1S(49:N1)*2/dy/dy;
  Q1bc(:,M1) = Q1bc(:,M1) + q1N/dy;
  Q1bc       = Q1bc(:);
  
  %...set up omega2
  q2    = zeros(N2,M2); q2  = q2(:);
  q2S =  zeros(N2,1);
  q2W =  zeros(1,M2);

%   q2N =  [zeros(N2-1,1);1];
%   q2E =  [zeros(1,M2-1),1];
 
  Q2bc       = zeros(N2,M2);
%   Q2bc(1,1:32)  = q2W(1:32)/dx;
%   Q2bc(1,33:N2)  = q2W(33:N2)*2/dx/dx;
%   Q2bc(N2,:) = q2E/dx;
%   Q2bc(:,1)  = Q2bc(:,1)  + q2S/dy;
%   Q2bc(:,M2) = Q2bc(:,M2) + q2N/dy;
  Q2bc       = Q2bc(:);
  
  f1 = ones(N1,M1);
%   f1(38:47,10:19) = 1;
  f1 = f1(:);
  
  f2 = ones(N2,M2);
  f2 = f2(:);
  
  
  sol = zeros(N,M);
  sol(1:48,1:32) = NaN;

  stopp = 1e-6;

  interfSp = zeros(16,1);
  interfEp = zeros(1,32);
  interfWp = zeros(1,32);
  
  % For plotting residuals
  Residuals = zeros(100,1);
  
  
count = 0;

q1old = q1;
q2old = q2;

%   while 1
errr = zeros(20,1);

for j=1:500
    count = count + 1;

    interfSpold = interfSp;
    interfEpold = interfEp;
    interfWpold = interfWp;

    q1    = PoissonSolve(q1old,f1-Q1bc,Rp1,Ap1,Pp1,N1*xLen/N,yLen/2);
    [Q2bc, interfWp] = setBoundary2p(q1,N1,M1,Q2bc,N2,M2);
    q2    = PoissonSolve(q2old,f2-Q2bc,Rp2,Ap2,Pp2,xLen/2,yLen);
    [Q1bc, interfSp, interfEp] = setBoundary1p(q2,N2,M2,Q1bc,N1,M1);

    GreenInterf = norm(interfWp-interfWpold);
    BlueInterf  = norm([(interfSp-interfSpold);(interfEp-interfEpold)']);
    Res = max([BlueInterf,GreenInterf]);
    errr(count) = get_error(q1,q2,N1,M1,N2,M2);
    % For plotting residuals
    Residuals(count) = Res;
%     Res

%     if Res < stopp
%         break;
%     end
end

  plot(1:length(errr),errr,'-o');
  
%   %...putting the sub-solutions together
% 
%   sol1 = reshape(q1,N1,M1);
%   sol2 = reshape(q2,N2,M2);
% 
%   %...graphics output
%   subplot(2,2,1)
%   sol(49:end,:)   = sol2;
%   sol(1:48,33:end) = sol1(1:48,:);  
%   soll2 = sol;
%   contourf(xx,yy,sol);
%   colorbar
%   title(['Contour-plot soution to' char(10) 'Poisson Equation with Omega 2 on top'])
%   axis([0 xLen 0 yLen])
%   axis equal
% 
%   subplot(2,2,3)
%   surf(xx,yy,sol);
%   title(['Surface-plot soution to' char(10) 'Poisson Equation with Omega 2 on top'])
%   axis tight
%   view(-35,35)
%     
%   subplot(2,2,2)
%   sol(49:end,:)   = sol2;
%   sol(1:N1,33:end) = sol1;
%   
%   contourf(xx,yy,sol);
%   colorbar
%   title(['Contour-plot soution to' char(10) 'Poisson Equation with Omega 1 on top'])
%   axis([0 xLen 0 yLen])
%   axis equal
%   
%   subplot(2,2,4)
%   surf(xx,yy,sol);
%   title(['Surface-plot soution to' char(10) 'Poisson Equation with Omega 1 on top'])
%   axis tight
%   view(-35,35)
% 
% %     figure(1)
% %   colormap(mmap)
% %   surf(xx,yy,u2D);
% %   shading interp;
% %   lighting phong;
% %   camlight headlight;
% %   material shiny;
% %   view(45,35)
% %   axis on
%   
%  % For plotting residuals
%  figure(2)
%  plot(1:count,Residuals(1:count),'-x');
%  title('Residuals in Poisson Equation')
%  ylabel('Residual')
%  xlabel('Iteration number in the while loop of Schwarz algorithm')
% 
%  
%  
%   figure(3)
%   contour(xx,yy,soll2-sol);
%   display(max(max(abs(soll2-sol))));