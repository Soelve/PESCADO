% This script determines spectrum of non-Hermitian Hamiltonian 
% (H with CAP). It also calculates the projections of the final state
% onty these eigenstates.

% Allocate
UmatComplex = zeros(N,N,lmax+1);
EnergyComplex = zeros(N,lmax+1);
Projections = zeros(N, lmax+1);

% Construct matrices with eigenvalues and states for non-Hermitian H
disp('Constructing the Non-Hermitian spectra.')
for L = 0:lmax
  % Hamiltonian for given L
  H_NH = T + diag(V) + L*(L+1)*diag(Sr) - 1i*diag(Gamma);
  % Eigenvalues and vectors
  [UmatNH EvalNH]=eig(H_NH);
  EvalNH=diag(EvalNH);
  % Normalize and invert
  UmatNH = UmatNH/sqrt(h);
  % NB: There is the option between LU decomposition ("inv")
  % and Moore Penrose pseudo inverse (pinv).
  UmatTilde = 1/h*(inv(UmatNH))';  % LU version
  %UmatTilde = 1/h*(pinv(UmatNH))';  % Pseudo inverse version    
  % Assign complex eigenenergies and eigenvectors to matrix
  EnergyComplex(:,L+1) = EvalNH;
  UmatComplex(:,:,L+1) = UmatNH;
  % Calculate projections
  Projections(:,L+1) = h * UmatTilde'*Psi(:,L+1);
  ProjectionsNoTilde(:,L+1) = h * UmatNH'*Psi(:,L+1);
  % Matrix elements between Herm. and NH eigenstates
  % BigMatWithGammaNH(:,:,L+1) = h * Umat(:,:,L+1)'*diag(Gamma)*UmatNH;
  % BigMatNoGammaNH(:,:,L+1) = h * Umat(:,:,L+1)'*UmatNH;
  % Write progress to screen
  %disp(['Done with L=',num2str(L),'.'])
end