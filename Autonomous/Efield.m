function F=Efield(E0,w,Ncycle,cep,t)

% Model for electric field in the dipole approximation:
% for 0 < t < Tpulse:
% E(t)=-A'(t) where 
% A(t) = E0/w sin^2(pi*t/Tpulse)*cos(wt + cep)
% Otherwise: E=0.
% E0: Maximum field strength
% w: Central angular frequency ("photon energy")
% Ncycle: Number of optical cycles
% cep: Carrier envelope phase
% t: time
% Alle quantities are given in atomic units.

Tpulse=Ncycle*2*pi/w;

%F=E0/w*(t>0).*(t<Tpulse).*...
%(2*pi/Tpulse*sin(pi*t/Tpulse).*cos(w*t+cep)-...
%w*(sin(pi*t/Tpulse)).^2.*sin(w*t+cep));

F=E0*(t>0).*(t<Tpulse).*...
(sin(pi*t/Tpulse).^2.*sin(w*t+cep)-...
pi/Tpulse/w*sin(2*pi*t/Tpulse).*cos(w*t+cep));