function [Rd,Pd] = getRPd(N,M)

%======================================
% set up a hierarchy of restriction
% and prolongation matrices (as cell
% arrays)
%======================================

  %...number of levels and initialization
  kk      = log(N)/log(2)-1;
  ll      = log(M)/log(2)-1;
  kl      = min(kk,ll);
  Rd      = cell(1,kl);
  Pd      = cell(1,kl);

  %...set up prolongation
  NN = N/2;
  MM = M/2;
  for i=1:kl
    PPx = sparse(2*NN,NN);
    for j=1:NN-1
      PPx(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
    end
    PPx(1,1)     = 0.5;
    PPx(2*NN,NN) = 0.5;
    NN = NN/2;

    PPy = sparse(2*MM,MM);
    for j=1:MM-1
      PPy(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
    end
    PPy(1,1)     = 0.5;
    PPy(2*MM,MM) = 0.5;
    MM = MM/2;

    Pd{i} = kron(PPy,PPx);
  end

  %...set up restriction (transpose of prolongation)
  for i=1:kl
    Rd{i} = transpose(Pd{i});
  end
