%--------------------------------------
% Driver routine for cavity flow
%--------------------------------------

  clear all
  close all

  addpath ~/Desktop/FOR_APS/NumPDE_Codes/colormaps/

  global Rd Ad Pd
  global Rp Ap Pp
  global xLen yLen
  global dx dy
  global N M
  global dt Re

  %...Reynolds number
  Re = 500;

  %...domain size
  xLen  = 1;
  yLen  = 1;

  %...resolution, fine-mesh spacing
  %   and number of cycles
  p     = 6;
  q     = 6;
  ncyc  = 10;

  N     = 2^p;
  M     = 2^q;

  x       = linspace(0,xLen,N+1);
  y       = linspace(0,yLen,M+1);
  xc      = (x(1:end-1)+x(2:end))/2;
  yc      = (y(1:end-1)+y(2:end))/2;
  [yy,xx] = meshgrid(yc,xc);
  dx      = xLen/N;
  dy      = yLen/M;
  dt      = min(dx,dy)/1.5;

  %...setup the hierarchy of R,A,P
  getRPd(N,M);
  getRPp(N,M);
  getAd(N,M);
  getAp(N,M);

  p0   = zeros(N,M); p0 = p0(:);
%   u1    = zeros(N,M); u1 = u1(:);
% %   v1    = zeros(N,M); v1 = u1(:);

  % EXTRA1
% %   v2 = u1;
%   u2 = u1;
  
  c1 = 1;
  c2 = 0;
  
  N1 = 48;
  M1 = 32;
  
  u1 = zeros(N1,M1); u1 = u1(:);
  c1u1 = c1*ones(size(u2));
  c2u1 = c2*ones(size(u2));

  u1S =  zeros(N1,1);
  u1N =  ones(N1,1);
  u1W =  zeros(1,M1);
  u1E =  zeros(1,M1);  
  
% %   v1S = zeros(N1,1);
% %   v1N = zeros(N1,1);
% %   v1W = zeros(1,M1);
% %   v1E = zeros(1,M1);
  
  N2 = 32; 
  M2 = 64; 
  
  u2 = zeros(N2,M2); u2 = u2(:);
  c1u2 = c1*ones(size(u2));
  c2u2 = c2*ones(size(u2));

  u2S =  zeros(N2,1);
  u2N =  ones(N2,1);
  u2W =  zeros(1,M2);
  u2E =  zeros(1,M2);

  x       = linspace(0,xLen,N+1);
  y       = linspace(0,yLen,M+1);
  xc      = (x(1:end-1)+x(2:end))/2;
  yc      = (y(1:end-1)+y(2:end))/2;
  [yy,xx] = meshgrid(yc,xc);
  dx      = xLen/N;
  dy      = yLen/M;
  dt      = min(dx,dy)/1.5;
  
% %   v2S = zeros(N2,1);
% %   v2N = zeros(N2,1);
% %   v2W = zeros(1,M2);
% %   v2E = zeros(1,M2);

%   Ubc      = zeros(N,M);
%   Ubc(1,:) = u1W*dt/2/Re/dx/dx;
%   Ubc(N,:) = u1E*dt/2/Re/dx/dx;
%   Ubc(:,1) = u1S*dt/2/Re/dy/dy;
%   Ubc(:,M) = u1N*dt/2/Re/dy/dy;
%   Ubc = Ubc(:);
% 
%   Vbc      = zeros(N,M);
%   Vbc(1,:) = v1W*dt/2/Re/dx/dx;
%   Vbc(N,:) = v1E*dt/2/Re/dx/dx;
%   Vbc(:,1) = vS*dt/2/Re/dy/dy;
%   Vbc(:,M) = v1N*dt/2/Re/dy/dy;
%   Vbc = Vbc(:);

  %extra
% %   sol = zeros(N1+N2,M2);
% %   sol(1:N1,1:M2-M1) = NaN; % for the plot
  
  for i=1:2000

    %...semi-Lagrangian advection
    u1    = SemiLagrAdvectNew(c1u1,c2u1,u1,u1S,u1N,u1W,u1E);
    temp = reshape(u1,N1,M1); % can obtain new u1E as u1(N1:N1:end)
    tempEast = temp(N1,:);
    u2W(M2-M1+1:end)= tempEast;%u1(end-M1+1:end);
    
    u2    = SemiLagrAdvectNew(c1u2,c2u2,u2,u2S,u2N,u2W,u2E);
    
    
%     v1    = SemiLagrAdvectNew(u1,v1,v1,vS,v1N,v1W,v1E);
    
%     v2W(1:M1)= v1(end- M1+1:end);
%     v2    = SemiLagrAdvectNew(u2,v2,v2,v2S,v2N,v2W,v2E);

%     sol(1:N2,:) = reshape(;
                    %     %...diffusion
                    %     unew    = DiffusionSolve(uast,Ubc);
                    %     vnew    = DiffusionSolve(vast,Vbc);
                    % 
                    %     %...computing divergence
                    %     [ux,uy] = Diff(unew,N,M,'D',uS,uN,uW,uE);
                    %     [vx,vy] = Diff(vnew,N,M,'D',vS,vN,vW,vE);
                    %     Div     = ux + vy; Div = Div(:); Div(1) = 0;
                    % 
                    %     %...solving for pressure
                    %     pnew    = PoissonSolve(p0,Div);
                    %     p0      = pnew;
                    % 
                    %     %...correcting velocities
                    %     [px,py] = Diff(pnew,N,M,'N',vS,vN,vW,vE);
                    %     u       = unew - px(:);
                    %     v       = vnew - py(:);

    if (mod(i,25)==0)
      k = 2;
      fprintf('time step %i \n',i)
      %...graphics output
%     [ux,uy] = Diff(u1,N,M,'D',u1S,u1N,u1W,u1E);
%     [vx,vy] = Diff(v1,N,M,'D',vS,v1N,v1W,v1E);
%     vort    = uy-vx;
%     vort2D  = reshape(vort,N,M);
%     u2D  = reshape(u1,N,M);
%     v2D  = reshape(v1,N,M);
%     vel  = sqrt(u2D.*u2D + v2D.*v2D);

      figure(1);
%     quiver(xx(1:k:end,1:k:end),yy(1:k:end,1:k:end),u2D(1:k:end,1:k:end),v2D(1:k:end,1:k:end),15/k)
%     contourf(xx,yy,vel,100,'EdgeColor','None')
      axis image
%     axis([0 xLen 0 yLen]);
      drawnow
    end
  end