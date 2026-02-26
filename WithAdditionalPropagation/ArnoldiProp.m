function PsiNew=LanczosProp(H0rad,Sr,Sl,rVector,...
Pz,AngularOp,AngularOp2,E0,w,Ncycle,cep,t,Psi,n,dt,h)

% This routine implements the Lanczos propagator
% and uses it to propagate the state Psi a time-step
% dt. n is the dimension of the Krylov subspace.
% The action of the Hamiltonian is provided by
% the subroutine HamMult.m.
% H0rad, Sr and Sl are related to the Hamiltonian.
% Psi is the state.
% dt is the time-step. (For imaginary time this should be
% set to -i times the actual dt in the call.)
% h is the radial increment in the grid.

% Allocation
[Rows Cols]=size(Psi);
Vbig=zeros(Rows,Cols,n);
Hmat = zeros(n+1,n);
T=zeros(n,n);

% Construct Arnoldi basis and Arnoldi Hamiltonian
InitialNorm=NormMatrixPsi(Psi,h); % In case of non-Hermcity: Perserver norm

Vbig(:,:,1)=Psi/InitialNorm;                            % First state
for k=2:n+1;
  V = Vbig(:,:,k-1);
  % Apply Hamiltonian
  U = HamMultUnpert(H0rad,Sr,Sl,V);
  % Velocity gauge
  U = U + HamMultPertVG(rVector,Pz,AngularOp,AngularOp2,...
      V,E0,w,Ncycle,cep,t);  
  % Length gauge
  %U = U + HamMultPertLG(rVector,Pz,AngularOp,AngularOp2,...
  %    V,E0,w,Ncycle,cep,t);  
  %
  % Remove component of V_{j-1}
  for j = 1:(k-1)
    Hmat(j,k-1)=InnerProduct(Vbig(:,:,j),U,h);
    U = U - Hmat(j,k-1)*Vbig(:,:,j);
  end
  Hmat(k,k-1) = NormMatrixPsi(U,h);
  % Stop iterating if singular
  if Hmat(k,k-1) > 1e-10
    Vbig(:,:,k) = U/Hmat(k,k-1);  
  else
    warning('The Krylov representation of H is near singular')
    break
  end  
end

% Exponentiate Hamiltonian to construct approximate propagator
Hmat = Hmat(1:n,:);
PsiArnoldi = expm(-1i*dt*Hmat);
% Apply to Psi - extract first column
PsiArnoldi = PsiArnoldi(:,1);

% Reconstruct propagated state
PsiNew=zeros(Rows,Cols);
for k=1:n;
  PsiNew = PsiNew + PsiArnoldi(k)*Vbig(:,:,k);
end

% Renormalize
PsiNew = PsiNew*InitialNorm;