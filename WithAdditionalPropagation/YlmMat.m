function Y=YlmMat(theta,L)

% The shpherical harmonic Y_{l,m}(\theta,\phi)
% with m=0. (It is \phi-independent.)

[a b]=size(theta);
Y=zeros(a,b);

Leg=legendre(L,cos(theta)); % Legendre polynomial

if L>0
  Y=Leg(1,:,:);
end

Y=sqrt((2*L+1)/(4*pi))*Y;
Y=reshape(Y,a,b);