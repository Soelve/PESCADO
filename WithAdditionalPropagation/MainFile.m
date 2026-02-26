% Initiate grid and operators (Inputs are given in a separate file)
clear

%Inputs                  % Provide input parameters (hard coded in separate script)
Inputs

SetUpH                  % Set up the Hamiltonian

% Construct initial state by imaginary time 
% (The inputs are hard coded within the subroutine)
disp('Constructing initial state by propagation in imaginary time')
ImaginaryTime
disp('Done')

% Propagate
disp('Propagating wave packet')
%PropagationLoop         % Call to the actual propagation
tic
PropagationLoopWithAnalysis
toc

% Present results
disp('Present results')
PresentResults
disp('Done with everything')