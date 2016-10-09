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
                
                %...Dirichlet b.c. on west
                PPx = sparse(2*NN,NN);
                for j=1:NN-1
                  PPx(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPx(1,1)     = 1;
                PPx(2*NN,NN) = 0.5;

                %...Neumann b.c.
                PPyN = sparse(2*MM,MM);
                for j=1:MM*3/4-1
                  PPyN(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPyN(1,1)     = 1;
                PPyN(2*MM,MM) = 1;

                %...Dirichlet b.c on south
                PPyDir = sparse(2*MM,MM);
                for j=MM*3/4:MM-1
                  PPyDir(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPyDir(1,1)     = 0.5;
                PPyDir(2*MM,MM) = 1;

                Pp{i} = kron(PPyN,PPx) + kron(PPyDir,PPx);

                NN = NN/2;
                MM = MM/2;

            end

  elseif (domain == 2)

            for i=1:kl
                PPx = sparse(2*NN,NN);
                for j=1:NN/2-1
                  PPx(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPx(1,1)     = 1;
                PPx(2*NN,NN) = 1;

                PPxDir = sparse(2*NN,NN);
                for j=NN/2:NN-1
                  PPxDir(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPxDir(1,1)     = 0.5;
                PPxDir(2*NN,NN) = 1;

                PPy = sparse(2*MM,MM);
                for j=1:MM-1
                  PPy(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
                end
                PPy(1,1)     = 1;
                PPy(2*MM,MM) = 1;

                Pp{i} =  kron(PPy,PPx)+kron(PPy,PPxDir);


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
