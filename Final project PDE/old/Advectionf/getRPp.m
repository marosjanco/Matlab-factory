function getRPp(N,M)
%======================================
% set up a hierarchy of restriction
% and prolongation matrices (as cell
% arrays) (Neumann version)
%======================================

  global Rp Ap Pp
  global ilevmin

  %...number of levels and initialization
  kk      = log(N)/log(2)-1;
  ll      = log(M)/log(2)-1;
  kl      = min(kk,ll);
  ilevmin = kl;
  Rp      = cell(1,kl);
  Pp      = cell(1,kl);

  %...set up prolongation
  NN = N/2;
  MM = M/2;
  for i=1:kl
    PPx = sparse(2*NN,NN);
    for j=1:NN-1
      PPx(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
    end
    PPx(1,1)     = 1;
    PPx(2*NN,NN) = 1;
    NN = NN/2;

    PPy = sparse(2*MM,MM);
    for j=1:MM-1
      PPy(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
    end
    PPy(1,1)     = 1;
    PPy(2*MM,MM) = 1;
    MM = MM/2;

    Pp{i} = kron(PPy,PPx);
  end

  %...set up restriction (transpose of prolongation)
  for i=1:kl
    Rp{i} = transpose(Pp{i});
  end
