%--------------------------------------
% Driver routine for cavity flow
%--------------------------------------

  clear all
  close all

  addpath ~/Desktop/FOR_APS/NumPDE_Codes/colormaps/

  global xLen yLen
  global dx dy
  global N M
  global dt Re

  %...Reynolds number
  Re = 500;


  %...resolution, fine-mesh spacing
  %   and number of cycles
  p     = 6;
  q     = 6;
  ncyc  = 10;

  c1 = 1;
  c2 = 0;
  
  N1 = 32;
  M1 = 32;
  
  u1 = zeros(N1,M1); u1 = u1(:);
  c1u1 = c1*ones(size(u1));
  c2u1 = c2*ones(size(u1));

  u1S =  zeros(N1,1);
  u1N =  ones(N1,1);
  u1W =  zeros(1,M1);
  u1E =  zeros(1,M1);  
  
  N2 = 32;
  M2 = 32;
  
  u2 = zeros(N2,M2); u2 = u2(:);
  c1u2 = c1*ones(size(u2));
  c2u2 = c2*ones(size(u2));

  u2S =  zeros(N2,1);
  u2N =  ones(N2,1);
  u2W =  zeros(1,M2);
  u2E =  zeros(1,M2);

  %...domain size
  xLen  = 1;
  yLen  = 0.5;
  
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
  
  f_init = @(x,y) exp( -100.*((x-0.3).^2+(y-0.2).^2));
  sol = f_init(xx,yy);

  u1 = sol(1:N1,:); u1 = u1(:);
  u2 = sol(N1+1:end,:); u2 = u2(:);
    
  utrue = sol; utrue = utrue(:);
  c1utrue = c1*ones(size(utrue));
  c2utrue = c2*ones(size(utrue));

  utrueS =  zeros(N,1);
  utrueN =  ones(N,1);
  utrueW =  zeros(1,M1);
  utrueE =  zeros(1,M1);

  for i=1:2000

    u1    = SemiLagrAdvect(c1u1,c2u1,u1,u1S,u1N,u1W,u1E,N1,M1);
    sol1 = reshape(u1,N1,M1);
    u2W= sol1(N1,:);
    
    u2W= (sol1(N1,:)+sol2(1,:))/2
    
    u2    = SemiLagrAdvect(c1u2,c2u2,u2,u2S,u2N,u2W,u2E,N2,M2);
    sol2 = reshape(u2,N2,M2);
    u1E = sol2(1,:);

    sol(1:N1,M2-M1+1:end) = sol1;
    sol(N1:end-1,:) = sol2;
    
    utrue    = SemiLagrAdvect(c1utrue,c2utrue,utrue,utrueS,utrueN,utrueW,utrueE,N,M);
    soltrue = reshape(utrue,N,M);
% %     when i have     sol(N1:end-1,:) = sol2;  0.0201
    compare = (sol(1:N-1,:) - soltrue(1:N-1,:));
    max(max(abs(compare)))

    if (mod(i,1)==0)
%       pause(0.05);
      k = 2;
      fprintf('time step %i \n',i)
      
      subplot(2,1,1)   
      %   colormap(mmap)
%       view(45,35)
      contourf(xx,yy,soltrue);
%       shading interp;
%       lighting phong;
      colorbar;
      camlight headlight;
%       material shiny;
      xlabel('x');
      ylabel('y');
      axis([0 xLen 0 yLen]);
      axis image;
      axis on;
      
      subplot(2,1,2)
      contourf(xx,yy,sol);
      colorbar;
      camlight headlight;
%       material shiny;
      xlabel('x');
      ylabel('y');
      axis([0 xLen 0 yLen]);
      axis image;
      axis on;
      drawnow
    end
  end