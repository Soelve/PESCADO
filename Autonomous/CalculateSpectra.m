% This script constructs box-noramlized basis of energy 
% eigen states for the Coulomb problem. 
% 

% Allocate
U=zeros(N,N,lmax+1);
Emat = zeros(N, lmax+1);

disp('Constructing the Hermitian spectrum.')
for L=0:lmax
  % Hamiltonian for given L
  H=T+diag(V) + L*(L+1)*diag(Sr);
  % Eigenvalues and vectors
  [Evect Eval]=eig(H);
  Eval=diag(Eval); [Eval, Eind]=sort(Eval);
  Evect=Evect(:,Eind);
  % Find lowest positive eigenenergy
  MinPos = Emat(min(find(Emat>0)));
  % Reassign minimal energy on grid
  if Emin < MinPos
      Emin = MinPos;
      disp(['New lower energy limit: ',num2str(Emin),...
          ', for l=',num2str(L)])
  end
  % Impose correct sign
  for n=1:N
    if Evect(2,n)<0
      Evect(:,n)=-Evect(:,n);
    end
  end
  % Normalize eivenvectors and assign to big array
  Umat(:,:,L+1)=Evect/sqrt(h);
  % Write eigenenergies to a matrix
  Emat(:,L+1)=Eval;
  disp(['Done with L=',num2str(L),'.'])
end