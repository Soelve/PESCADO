function State=HamMultPert(rVector,Pz,AngularOp,...
    AngularOp2,Psi,E0,w,Ncycle,cep,t)

% Function which calculates the action of the perturbation 
% on the wave functon. The wave function is given as a matrix
% in which the rows correspond to the radial degrees of freedom, 
% and the columns correspond to the angular channel in the
% dipole approximation.

Aux = 0*Psi;
State = Psi*AngularOp;
State = Pz*State;
Aux = Psi*AngularOp2;
Aux = -1i*diag(1./rVector)*Aux;
State = State + Aux;
%State = State - Aux;                 % Feil????
State= Afield(E0,w,Ncycle,cep,t)*State;
