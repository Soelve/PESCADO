function F=AfieldInt(E0,w,Ncycle,cep,t)

% Model for electric field in the dipole approximation:
% for 0 < t < Tpulse:
% A(t)=E0/w * sin^2(pi*t/Tpulse)*cos(wt + cep)
% Otherwise: A=0.
%
% This function calculates the integral 
% int_0^t dt' A(t').
%
% E0: Maximum field strength
% w: Central angular frequency ("photon energy")
% Ncycle: Number of optical cycles
% cep: Carrier envelope phase
% t: time
% Alle quantities are given in atomic units.

Tpulse=Ncycle*2*pi/w;

F=E0/w/4*(t>0).*(t<Tpulse).*(...
Tpulse*sin(2*pi*t/Tpulse-w*t-cep)/(Tpulse*w-2*pi)-...
Tpulse*sin(2*pi*t/Tpulse+w*t+cep)/(Tpulse*w+2*pi)+...
2*cos(cep)*sin(w*t)/w+2*sin(cep)*cos(w*t)/w);

F=F-E0/w/4*(t>0).*(t<Tpulse).*(...
Tpulse*sin(-cep)/(Tpulse*w-2*pi)-...
Tpulse*sin(+cep)/(Tpulse*w+2*pi)+...
2*sin(cep)/w);