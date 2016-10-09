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
  ncyc  = 10;

  



  u    = zeros(N,M); u = u(:);

   %extra
   
  u1 = u;
  u2 = u;
  
  N1 = 48;
  M1 = 32;
  
  u1S =  zeros(N1,1);
  u1N =  ones(N1,1);
  u1W =  zeros(1,M1);
  u1E =  zeros(1,M1);

  N2 = 32;
  M2 = 64;
  
  u2S =  zeros(N2,1);
  u2N =  ones(N2,1);
  u2W =  zeros(1,M2);
  u2E =  zeros(1,M2);

  Len = N1+N2;
  x       = linspace(0,Len,Len+1);
  y       = linspace(0,Len,Len+1);
  xc      = (x(1:end-1)+x(2:end))/2;
  yc      = (y(1:end-1)+y(2:end))/2;
  [yy,xx] = meshgrid(yc,xc);
  dx      = xLen/N;
  dy      = yLen/M;
  dt      = min(dx,dy)/1.5;
  

  
  %extra
  sol = zeros(N1+N2,M2);
  sol(N2+1:end,M1:end) = NaN;
  
  for i=1:2000

    %...semi-Lagrangian advection
    u1    = SemiLagrAdvectNew(zeros(N1*M1),ones(N1*M1),u1,u1S,u1N,u1W,u1E);
    uu1 = reshape(u1,[N1,M1]);
    u2W(1:M1) = uu1(end,:);
    u2    = SemiLagrAdvectNew(zeros(N2*M2),ones(N2*M2),u2,u2S,u2N,u2W,u2E);

    sol(1:N2,:) = reshape(u2,[N2,M2]);
    sol(N2 +1:end,M1) = reshape(u1,[N1,M1]);
    
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

    if (mod(i,5)==0)
      k = 2;
      fprintf('time step %i \n',i)
      %...graphics output
      co
      
      
      
      
      
%       [ux,uy] = Diff(u1,N,M,'D',u1S,u1N,u1W,u1E);
%       [vx,vy] = Diff(v1,N,M,'D',vS,v1N,v1W,v1E);
%       vort    = uy-vx;
%       vort2D  = reshape(vort,N,M);
%       u2D  = reshape(u1,N,M);
%       v2D  = reshape(v1,N,M);
%       vel  = sqrt(u2D.*u2D + v2D.*v2D);
      
        %...graphics output

      mmap = cbrewer('div','RdBu',64);
      figure(1)
      colormap(mmap)
      surf(xx,yy,sol);
      shading interp;
      lighting phong;
      camlight headlight;
      material shiny;
      view(45,35)
  
%       figure(1);
% %       quiver(xx(1:k:end,1:k:end),yy(1:k:end,1:k:end),u2D(1:k:end,1:k:end),v2D(1:k:end,1:k:end),15/k)
%      contourf(xx,yy,vel,100,'EdgeColor','None')
%       axis image
%       axis([0 xLen 0 yLen]);
      drawnow
    end
  end
