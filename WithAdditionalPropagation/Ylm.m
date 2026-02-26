function Y=YlmMat(theta,L)

% The shpherical harmonic Y_{l,m}(\theta,\phi)
% with m=0. (It is \phi-independent.)

Leg=legendre(L,cos(theta)); % Legendre polynomial

%Leg=Leg(L+1,:);
Leg=Leg(1,:);

Y=sqrt((2*L+1)/(4*pi))*Leg;