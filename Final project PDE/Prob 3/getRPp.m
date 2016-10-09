function [Rp,Pp] = getRPp(N,M,domain)

%======================================
% set up a hierarchy of restriction
% and prolongation matrices (as cell
% arrays) (Neumann version)
%======================================

  %...number of levels and initialization
  kk      = log(N)/log(2)-1;
  ll      = log(M)/log(2)-1;
  kl      = min(kk,ll);
  Rp      = cell(1,kl);
  Pp      = cell(1,kl);

  %...set up prolongation
  NN = N/2;
  MM = M/2;
  
   if (domain == 1)
       
            for i=1:kl
                PPx = sparse(2*NN,NN);

                for j=1:NN-1
                  PPx(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPx(1,1)     = 1;
                PPx(2*NN,NN) = 0.5;
               

                PPy = sparse(2*MM,MM);
                for j=1:MM-1
                  PPy(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPy(1,1)     = 1;
                PPy(2*MM,MM) = 1;
                
                Ptemp = kron(PPy,PPx);
                Ptemp(1:2*NN,NN*3/4+1:NN) = Ptemp(1:2*NN,NN*3/4+1:NN)/2;
                Pp{i} = Ptemp;
                
                NN = NN/2;
                MM = MM/2;
            end

  elseif (domain == 2)

            for i=1:kl
                PPx = sparse(2*NN,NN);
                for j=1:NN-1
                  PPx(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPx(1,1)     = 1;
                PPx(2*NN,NN) = 1;

                PPy = sparse(2*MM,MM);
                for j=1:MM-1
                  PPy(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPy(1,1)     = 1;
                PPy(2*MM,MM) = 1;

                PPxDir = PPx;
                PPxDir(1,1) = 0.5;

                PPySouth = zeros(2*MM,MM);
                PPySouth(1:MM/2,:) = PPy(1:MM/2,:);
                PPyNorth = PPy - PPySouth;
                
                Pp{i} = kron(PPySouth,PPx) + kron(PPyNorth,PPxDir);
    
                NN = NN/2;
                MM = MM/2;

            end

  else
     fprintf('\nError:\nIncorrect parameter for domain label.\n');
     return; 
   end
      
  %...set up restriction (transpose of prolongation)
  for i=1:kl
    Rp{i} = transpose(Pp{i});
  end
