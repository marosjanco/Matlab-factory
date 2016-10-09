%--------------------------------------
% Driver routine for cavity flow
%--------------------------------------
% 
%   clear all
%   close all

  global xLen yLen
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

  
  
  %...Dimensions of omega_1
      N1    = 64;  %Diffusion/Poisson
      N1adv = 48;  %Advection
      M1    = 32;  %Advection/Diffusion/Poisson
            
  %...Dimensions of omega_2 for Advection/Diffusion/Poisson
      N2 = 32;
      M2 = 64;
 
  %...setup 2D velocity fields
      u2D = [NaN*zeros(N1adv,M-M1),zeros(N1adv,M1);zeros(N2,M2)];
      v2D = u2D;
  
  %...Initializiation for Advection         
      u1_2D = zeros(N1adv,M1);
      v1_2D = zeros(N1adv,M1);

      u2_2D = zeros(N2,M2);
      v2_2D = zeros(N2,M2);
  
  %...Dirichlet external boundaries
      %...for u
          %...domain 1
              %Diffusion/Poisson
              u1S =  zeros(N1,1);
              u1N =  ones(N1,1);
              u1W =  zeros(1,M1);
              u1E =  zeros(1,M1);

              %Advection
              u1Sadv =  zeros(N1adv,1);
              u1Nadv =  ones(N1adv,1);
              u1Wadv =  zeros(1,M1);
              u1Eadv =  zeros(1,M1);

          %...domain 2 
              %Diffusion/Poisson
              u2S =  zeros(N2,1);
              u2N =  ones(N2,1);
              u2W =  zeros(1,M2);
              u2E =  zeros(1,M2);

              %Advection
              u2Sadv =  zeros(N2,1);
              u2Nadv =  ones(N2,1);
              u2Wadv =  zeros(1,M2);
              u2Eadv =  zeros(1,M2);

      %...for v   
          %...domain 1
              %Diffusion/Poisson
              v1S =  zeros(N1,1);
              v1N =  zeros(N1,1);
              v1W =  zeros(1,M1);
              v1E =  zeros(1,M1);

              %Advection
              v1Sadv =  zeros(N1adv,1);
              v1Nadv =  zeros(N1adv,1);
              v1Wadv =  zeros(1,M1);
              v1Eadv =  zeros(1,M1);
              
          %...domain 2  
              %Diffusion/Poisson
              v2S = zeros(N2,1);
              v2N = zeros(N2,1);
              v2W = zeros(1,M2);
              v2E = zeros(1,M2); 

              %Advection
              v2Sadv =  zeros(N2,1);
              v2Nadv =  zeros(N2,1);
              v2Wadv =  zeros(1,M2);
              v2Eadv =  zeros(1,M2);

  %...Embedding the Dirichlet boundaries
  %...... for u        
      U1bc       = zeros(N1,M1);
      U1bc(1,:)  = u1W*2*dt/Re/dx/dx;
      U1bc(N1,:) = u1E*2*dt/Re/dx/dx;
      U1bc(:,1)  = U1bc(:,1)  + u1S*2*dt/Re/dy/dy;
      U1bc(:,M1) = U1bc(:,M1) + u1N*2*dt/Re/dy/dy;
      U1bc       = U1bc(:);

      U2bc       = zeros(N2,M2);
      U2bc(1,:)  = u2W*2*dt/Re/dx/dx;
      U2bc(N2,:) = u2E*2*dt/Re/dx/dx;
      U2bc(:,1)  = U2bc(:,1)  + u2S*2*dt/Re/dy/dy;
      U2bc(:,M2) = U2bc(:,M2) + u2N*2*dt/Re/dy/dy;
      U2bc       = U2bc(:);

  %...... for v
      V1bc       = zeros(N1,M1);
      V1bc(1,:)  = v1W*2*dt/Re/dx/dx;
      V1bc(N1,:) = v1E*2*dt/Re/dx/dx;
      V1bc(:,1)  = V1bc(:,1)  + v1S*2*dt/Re/dy/dy;
      V1bc(:,M1) = V1bc(:,M1) + v1N*2*dt/Re/dy/dy;
      V1bc       = V1bc(:);

      V2bc       = zeros(N2,M2);
      V2bc(1,:)  = v2W*2*dt/Re/dx/dx;
      V2bc(N2,:) = v2E*2*dt/Re/dx/dx;
      V2bc(:,1)  = V2bc(:,1)  + v2S*2*dt/Re/dy/dy;
      V2bc(:,M2) = V2bc(:,M2) + v2N*2*dt/Re/dy/dy;
      V2bc       = V2bc(:);
  
  %...Neumann external boundaries (for velocity corrections)
      S1_Neum =  zeros(N1adv,1);
      N1_Neum =  zeros(N1adv,1);
      W1_Neum =  zeros(1,M1);
      E1_Neum =  zeros(1,M1);

      S2_Neum =  zeros(N2,1);
      N2_Neum =  zeros(N2,1);
      W2_Neum =  zeros(1,M2);
      E2_Neum =  zeros(0,M2);

  %...setup the pressure
      p10   = zeros(N1*M1,1);
      p20   = zeros(N2*M2,1);

  %...setup for Poisson boundaries
  P1bc = zeros(N1*M1,1);
  P2bc = zeros(N2*M2,1);
  
  %...Initialize the green and blue interfaces (blue is furthur splitted)
      interfS = zeros(16,1);
      interfE = zeros(1,32);
      interfW = zeros(1,32);
      
  %...setup the hierarchy of R,A,P for Diffusion
  %...... for omega_1
      [Rd1,Pd1]  = getRPd(N1,M1);
      [AAA1,Ad1] = getAd(N1,M1,Rd1,Pd1);
  %...... for omega_2
      [Rd2,Pd2]  = getRPd(N2,M2);
      [AAA2,Ad2] = getAd(N2,M2,Rd2,Pd2);
  
  %...setup the hierarchy of R,A,P for Poisson
  %...... for omega_1
      [Rp1,Pp1] = getRP2p(N1,M1,1);
      Ap1       = getAp(N1,M1,1,Rp1,Pp1);
  %...... for omega_2
      [Rp2,Pp2] = getRP2p(N2,M2,2);
      Ap2       = getAp(N2,M2,2,Rp2,Pp2);  
  
  %...setup stopping criteria
  stopd = 1e-7; %for Diffusion
  stopp = 1e-7; %for Poisson

  u1_2D = zeros(N1adv,M1);
  v1_2D = zeros(N1adv,M1);

  u2_2D = zeros(N2,M2);
  v2_2D = zeros(N2,M2);
 
  for i=1:2000

    %...semi-Lagrangian advection
    
    %...... extracting subsolutions for advection
    
        u1adv = u1_2D(:);
        v1adv = v1_2D(:);
        u2adv = u2_2D(:);
        v2adv = v2_2D(:);
        
    %...... advecting u
        
        %...Updating the boundaries between omega_1 and omega_2
        u1Eadv  = (u1_2D(N1adv,:)+u2_2D(1,33:end))/2;
        u2Wadv(33:end) = u1Eadv;

        %...Advecting omega_1 and omega_2 separately
        u1ast  = SemiLagrAdvect(u1adv,v1adv,u1adv,u1Sadv,u1Nadv,u1Wadv,u1Eadv,N1adv,M1);   
        u2ast   = SemiLagrAdvect(u2adv,v2adv,u2adv,u2Sadv,u2Nadv,u2Wadv,u2Eadv,N2,M2);

    %...... advecting v
    
        %...Updating the boundaries between omega_1 and omega_2
        v1Eadv  = (v1_2D(N1adv,:)+v2_2D(1,33:end))/2;
        v2Wadv(33:end) = v1Eadv;
        
        %...Advecting omega_1 and omega_2 separately
        v1ast  = SemiLagrAdvect(u1adv,v1adv,v1adv,v1Sadv,v1Nadv,v1Wadv,v1Eadv,N1adv,M1);   
        v2ast   = SemiLagrAdvect(u2adv,v2adv,v2adv,v2Sadv,v2Nadv,v2Wadv,v2Eadv,N2,M2);

        %...... extracting subsolutions for Diffusion
    
        u2_2D = reshape(u2ast,N2,M2);
        u1_2D = [reshape(u1ast,N1adv,M1);u2_2D(1:16,33:end)];
        v2_2D = reshape(v2ast,N2,M2);
        v1_2D = [reshape(v1ast,N1adv,M1);v2_2D(1:16,33:end)];
        
        u1ast = u1_2D(:);
        v1ast = v1_2D(:);

    %...diffusion
    %......diffusing u

          while 1      
            interfSold = interfS;
            interfEold = interfE;
            interfWold = interfW;

            u1new    = DiffusionSolve(u1ast,U1bc,AAA1,Rd1,Ad1,Pd1);
            [U2bc, interfW] = setBoundary2d(u1new,N1,M1,U2bc,N2,M2,u2N(1));
            u2new    = DiffusionSolve(u2ast,U2bc,AAA2,Rd2,Ad2,Pd2);
            [U1bc, interfS, interfE] = setBoundary1d(u2new,N2,M2,U1bc,N1,M1,u1N(end));

            NormGreen = norm(interfW-interfWold);
            NormBlue  = norm([(interfS-interfSold);(interfE-interfEold)']);
            Res = max([NormBlue,NormGreen]);
           %this nees to be the same!!!!
%           pp1 = reshape(u1new,N1,M1);
%           pp2 = reshape(u2new,N2,M2);
%           mx1 = pp1(49:end,:);
%           mx2 = pp2(1:16,33:end);
%           subplot(1,2,1);
%           contourf(mx1);
%           subplot(1,2,2);
%           contourf(mx2);pause(3);
            if (Res < stopd)
                break;
            end
          end
          
%         soltemp2 = reshape(u2new,N2,M2);
%         soltemp1 = reshape(u1new,N1,M1);
%       
%         soltemp1(57:end,:) = soltemp2(9:16,33:end);
%         soltemp2(1:8,33:end) = soltemp1(49:56,:);
%       
%         u1new = soltemp1(:);
%         u2new = soltemp2(:);
% 
    %......diffusing v

          while 1      
            interfSold = interfS;
            interfEold = interfE;
            interfWold = interfW;

            v1new    = DiffusionSolve(v1ast,V1bc,AAA1,Rd1,Ad1,Pd1);
            [V2bc, interfW] = setBoundary2d(v1new,N1,M1,V2bc,N2,M2,v2N(1));
            v2new    = DiffusionSolve(v2ast,V2bc,AAA2,Rd2,Ad2,Pd2);
            [V1bc, interfS, interfE] = setBoundary1d(v2new,N2,M2,V1bc,N1,M1,v1N(end));

            NormGreen = norm(interfW-interfWold);
            NormBlue  = norm([(interfS-interfSold);(interfE-interfEold)']);
            Res = max([NormBlue,NormGreen]);
            

            if (Res < stopd)
                break;
            end
          end
%           
%         soltemp2 = reshape(v2new,N2,M2);
%         soltemp1 = reshape(v1new,N1,M1);
%       
%         soltemp1(57:end,:) = soltemp2(9:16,33:end);
%         soltemp2(1:8,33:end) = soltemp1(49:56,:);
%       
%         v1new = soltemp1(:);
%         v2new = soltemp2(:);

    %...computing divergence
    
      u1_2D  = reshape(u1new,N1,M1); 
      v1_2D  = reshape(v1new,N1,M1); 
      
%       u2_2D  = reshape(u2new,N2,M2);
%       v2_2D  = reshape(v2new,N2,M2);
      
            
    %... divergence for omega_1
          uu1E = (u1_2D(N1adv,:)+u1_2D(N1adv+1,:))/2;
          vv1E = (v1_2D(N1adv,:)+v1_2D(N1adv+1,:))/2;
          
          u1_2D = u1_2D(1:N1adv,:); u1new = u1_2D(:);
          v1_2D = v1_2D(1:N1adv,:); v1new = v1_2D(:);
          
          [ux,~] = Diff(u1new,N1adv,M1,'D',u1Sadv,u1Nadv,u1Wadv,uu1E,1);
          [~,vy] = Diff(v1new,N1adv,M1,'D',v1Sadv,v1Nadv,v1Wadv,vv1E,1);

          Div1     = ux + vy; Div1 = Div1(:); Div1(1) = 0;

    %... divergence for omega_2
          uu2W = [u2Wadv(1:32),uu1E];
          vv2W = [v2Wadv(1:32),vv1E];
          
          [ux,~] = Diff(u2new,N2,M2,'D',u2Sadv,u2Nadv,uu2W,u2Eadv,2);
          [~,vy] = Diff(v2new,N2,M2,'D',v2Sadv,v2Nadv,vv2W,v2Eadv,2);
          Div2     = ux + vy; Div2 = Div2(:); Div2(1) = 0;

    %...solving for pressure
    
          P1 = reshape(p10,N1,M1);
          P2 = reshape(p20,N2,M2);
          
          P1bc = reshape(P1bc,N1,M1);
          P2bc = reshape(P2bc,N2,M2);
          
          P1bc(N1adv+1:end,1) = (P2(1:16,32)+P2(1:16,33))/2;
          P1bc(end,:)    = (P2(16,33:end)+P2(17,33:end))/2;
          
          P2bc(1,33:end) = (P1(N1adv,:)+P1(N1adv+1,:))/2;
          
          P1bc = P1bc(:);
          P2bc = P2bc(:);
          
          Div1 = reshape(Div1,N1adv,M1);
          Div2 = reshape(Div2,N2,M2);
          Div = [NaN*zeros(N1adv,M-M1),Div1;Div2];
          Div1 = Div(1:N1,33:end);
          Div1 = Div1(:); Div2 = Div2(:);

          
          c = 0;
          while 1
            
            c = c+1;
            interfSold = interfS;
            interfEold = interfE;
            interfWold = interfW;

            p1new    = PoissonSolve(p10,Div1-P1bc,Rp1,Ap1,Pp1);
            [P2bc, interfW] = setBoundary2p(p1new,N1,M1,P2bc,N2,M2);
            p2new    = PoissonSolve(p20,Div2-P2bc,Rp2,Ap2,Pp2);
            [P1bc, interfS, interfE] = setBoundary1p(p2new,N2,M2,P1bc,N1,M1);
            
            NormGreen = norm(interfW-interfWold);
            NormBlue  = norm([(interfS-interfSold);(interfE-interfEold)']);
            Res = max([NormBlue,NormGreen]);
            
%           this nees to be the same ! ! ! !
%           pp1 = reshape(p1new,N1,M1);
%           pp2 = reshape(p2new,N2,M2);
%           mx1 = pp1(49:end,:);
%           mx2 = pp2(1:16,33:end);
%           subplot(1,2,1);
%           contourf(mx1);
%           subplot(1,2,2);
%           contourf(mx2);pause(0.2);
%             display(c);
%             display(Res);
%                       display(max(max(abs(mx1-mx2))));

            if Res < stopp
                break;
            end
          end
          
          p10 = p1new;
          p20 = p2new;

          p1_2D = reshape(p1new,N1,M1);
          p2_2D = reshape(p2new,N2,M2);
                    
    %...correcting velocities
    %...... correcting omega_1 velocities
     
          E1_Neum = (p1_2D(N1adv,:)+p1_2D(N1adv+1,:))/2;
          
          pressure = [NaN*zeros(N1adv,M-M1),p1_2D(1:N1adv,:);p2_2D];
          p1_2D = pressure(1:N1adv,33:end);
          p1new = p1_2D(:);
          
      [px,py] = Diff(p1new,N1adv,M1,'N',S1_Neum,N1_Neum,W1_Neum,E1_Neum,1);
      u1       = u1new - px(:);
      v1       = v1new - py(:);

    %...... correcting omega_2 velocities        
      W2_Neum(33:end) = E1_Neum;
      
      [px,py] = Diff(p2new,N2,M2,'N',S2_Neum,N2_Neum,W2_Neum,E2_Neum,2);
      u2      = u2new - px(:);
      v2      = v2new - py(:);

      u1_2D  = reshape(u1,N1adv,M1);
      v1_2D  = reshape(v1,N1adv,M1);
      
      u2_2D  = reshape(u2,N2,M2);
      v2_2D  = reshape(v2,N2,M2);

%     if (mod(i,20)==0)
            if (i==20)

      k = 2;
      fprintf('time step %i \n',i)
      %...graphics output
            
      u2D= [NaN*zeros(N1adv,M-M1),u1_2D;u2_2D];
      v2D= [NaN*zeros(N1adv,M-M1),v1_2D;v2_2D];

      figure(1);
      quiver(xx(1:k:end,1:k:end),yy(1:k:end,1:k:end),u2D(1:k:end,1:k:end),v2D(1:k:end,1:k:end),12/k)
      title(sprintf('Solution at time-step = %i',i))
      axis image
      axis([0 xLen 0 yLen]);
      drawnow
      break
    end
  end
  
  
