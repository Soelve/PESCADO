function dPdO=AngularDist(theta,ThetaMatrix,Lmax)

% Function which calculates dP/d\Omega
% from an Hermitian matrix whith all the 
% weights i the incoherent sum of
% Spherical Harmonics-producs

% Initiate
dPdO=0;

% Loop over the matrix of weights
for L=0:Lmax
  for LL=0:Lmax
    NewTerm=ThetaMatrix(L+1,LL+1)*Ylm(theta,L).'.*...
        conj(Ylm(theta,LL)).';
    dPdO=dPdO+NewTerm;
  end
end

% Ensure reality
dPdO=real(dPdO);