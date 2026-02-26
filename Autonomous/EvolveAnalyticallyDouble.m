% This script solve the TDSE for t>T, after the laser
% pulse, when the Hamiltonian is time-independent.
% (It is also non-Hermitian).

% Allocate
BigUmat = zeros(N, N, lmax+1);
BigCmplxEmat = zeros(N, lmax+1);
Cmat = zeros(N, lmax+1);

% Loop over L and find eivenvalues
disp('Finding the eigenstates of the non-Hermitian Hamiltonian:')
for L = 0:lmax
  % Write L to screen
  disp(['L=',num2str(L)])
  % Non-Hermitian Hamiltonian for given l 
  Ham = H0rad + L*(L+1)*diag(Sr) -1i*diag(Gamma);
  % Eivenectors ana values
  [UmatNH CmplxE] = eig(Ham);
  % Sort and normalize
  CmplxE = diag(CmplxE); [CmplxE Inds] = sort(CmplxE); 
  UmatNH = UmatNH(:, Inds); UmatNH = UmatNH/sqrt(h);
  % Store
  BigUmat(:,:, L+1) = UmatNH;
  BigCmplxEmat(:, L+1) = CmplxE.';
  % Projections
  UmatTilde = 1/h*(inv(UmatNH))';
  Cmat(:,L+1) = h*UmatTilde'*PsiT(:,L+1);
end

% Cutoff in imaginary part
MaxImagE = -1e-5;

% Set up reminder of zeta
disp('Resolve the dynamics from t=T to infinity:')
DensMatDoubleAfter = zeros(N, N, lmax+1, lmax+1);
for L = 0:lmax
  for Lp = 0:lmax
    % Write L and Lp to screen
    disp(['L=',num2str(L),', Lp=',num2str(Lp)])
    for nL = 1:N
      En = BigCmplxEmat(nL,L+1);
      if imag(En) < MaxImagE
        Phi_n = BigUmat(:,nL,L+1);
        Phi_n = Gamma.*Phi_n;
        Cn = Cmat(nL, L+1);
        for mLp = 1:N
          Em = BigCmplxEmat(mLp, Lp+1);
          Cm = Cmat(mLp, Lp+1);
          Phi_m = BigUmat(:,mLp,Lp+1);
          DensMatDoubleAfter(:,:, L+1, Lp+1) = ...
          DensMatDoubleAfter(:,:,L+1,Lp+1) + ...
          1/(En-conj(Em))*Cn*conj(Cm)*Phi_n*Phi_m';
        end
      end
    end
  end
end
% Multiply -i and divide by h^2 - for proper normalization
DensMatDoubleAfter = -1i*DensMatDoubleAfter;