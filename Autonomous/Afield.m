function F=Afield(E0,w,Ncycle,cep,t)

% Model for electric field in the dipole approximation:
% for 0 < t < Tpulse:
% E(t)=E0 sin^2(pi*t/Tpulse)*cos(wt + cep)
% Otherwise:    E=0
% E0:           Maximum field strength
% w:            Central angular frequency ("photon energy")
% Ncycle:       Number of optical cycles
% cep:          Carrier envelope phase
% t:            time
% Alle quantities are given in atomic units.

% Single pulse
Tpulse=Ncycle*2*pi/w;
F=E0/w*(t>0).*(t<Tpulse).*(sin(pi*t/Tpulse)).^2.*cos(w*t+cep);

% Two pulses with a lag
%tau=5;
%F1=E0/w*(t>0).*(t<Tpulse).*(sin(pi*t/Tpulse)).^2.*cos(w*t+cep);
%t2=t-Tpulse-tau;
%F2=E0/w*(t2>0).*(t2<Tpulse).*(sin(pi*t2/Tpulse)).^2.*sin(w*t2+cep);
%
%F=F1+F2;