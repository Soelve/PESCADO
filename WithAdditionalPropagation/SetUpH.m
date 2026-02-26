% Skript som set opp H0

%
% Kinetic energy
%
rVector = linspace(0,Rmax,N+2)';
rVector = rVector(2:(N+1));
h = rVector(2)-rVector(1); 
T=zeros(N,N);
% Asymmetric five point formula
T(1,1:4)=[-20 6 4 -1];
% The rest: Symmetric five point formula (truncated in the ends)
T(2,1:4)=[16 -30 16 -1];
T(N-1,(N-3):N)=[-1 16 -30 16];
T(N,(N-2):N)=[-1 16 -30];
for nn=3:(N-2)
  T(nn,(nn-2):(nn+2))=[-1 16 -30 16 -1];
end
T=-1/2*T/12/h^2;

% Coulomb potential, V=-1/r
%V=-1./rVector;
% Yukawa potential
V=-exp(-Alpha*rVector)./rVector;

% Absorber
Gamma=eta*(rVector>Onset).*(rVector-Onset).^AbsPower;

% Total Hamiltonian - without the sentrifugal term 
H0rad=T+diag(V);
%clear T
% Sparsify the unperturbed, radial Hamiltonian
H0rad=sparse(H0rad);

% Sentrifugal term; +l(l+1)/2r^2
Sr=1/2./rVector.^2;
Lvector=0:lmax;
Sl=Lvector.*(Lvector+1);

% With Psi being an N \times (lmax+1) matrix:
% H \Psi is represented as
%(T + diag(V)) Psi + diag(Sr) Psi diag(Sl)

%
% Perturbation: E(t) diag(rVector)*Psi*AngularOp
%
AngularVector=(Lvector+1)./sqrt((2*Lvector+1).*(2*Lvector+3));
%AngularVector=(Lvector)./sqrt((2*Lvector-1).*(2*Lvector+1)); % FEIL!!!
AngularVector=AngularVector(1:(end-1));
AngularOp = diag(AngularVector,1); AngularOp=AngularOp+AngularOp.';
AngularOp2= diag(-(Lvector(1:(end-1))+1).*AngularVector,1); 
AngularOp2=AngularOp2-AngularOp2.';
%AngularOp=sparse(AngularOp);
%AngularOp2=sparse(AngularOp2);

% For velocity gauge: A(t) p_z*Psi*AngularOp 
% where p_z is simply -i d/dr \cos \theta 
% (in the dipole approximation)
%
% Three point formula:
%Pz=diag(ones(1,N-1),1); Pz=Pz-Pz.'; Pz=-1i*Pz/2/h;
%Pz=sparse(Pz);
% Five point formula:
Pz=zeros(N,N);              % Allocate
Pz(1,1:4)=[-10 18 -6 1];    % Asymmetric formula (Psi(0)=0)
% The rest: Symmetric five point formula
Pz(2,1:4)=[-8 0 8 -1];
Pz(N-1,(N-3):N)=[1 -8 0 8];
Pz(N,(N-2):N)=[1 -8 0];
for n=3:(N-3)
  Pz(n,(n-2):(n+2))=[1 -8 0 8 -1];
end
Pz=-1i*Pz/12/h;             % Introduce proper factors
Pz=sparse(Pz);              % Sparsify