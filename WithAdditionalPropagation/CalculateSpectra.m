% This script constructs box-noramlized basis of energy 
% eigen states for the Coulomb problem. 
% 

% Allocate
U=zeros(N,N,lmax+1);

disp('Constructing the spectrum.')
for L=0:lmax
  H=T+diag(V) + L*(L+1)*diag(Sr);
  [Evect Eval]=eig(H);
  Eval=diag(Eval); [Eval, Eind]=sort(Eval);
  Evect=Evect(:,Eind);
  for n=1:N
    if Evect(2,n)<0
      Evect(:,n)=-Evect(:,n);
    end
  end
  Umat(:,:,L+1)=Evect/sqrt(h);
  Emat(:,L+1)=Eval;
  disp(['Done with l= ',num2str(L),'.'])
end
%clear T Evect Eval

