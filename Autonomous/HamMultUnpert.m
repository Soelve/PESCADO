function State=HamMultUnpert(H0rad,Sr,Sl,Psi)

% Function which calculates the action of the
% Hamiltonian on the state matrix Psi.

% With absorber
%State= H0rad*Psi + diag(Sr)*Psi*diag(Sl)-1i*diag(Gamma)*Psi;
% Without absorber
State= H0rad*Psi + diag(Sr)*Psi*diag(Sl);