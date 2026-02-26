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
% Timing
t0 = cputime;
PropagationLoopWithAnalysis
CalulationTime = cputime - t0
%return

% Determine spectrum of non-Hermitian Hamiltonian - and load the analytical
% Coulomb waves
if EnergySpectrum | DoublyDiffSpectrum | AngularSpectrum
  CalculateSpectra
  CalculateNonHermSpectra
end    

% Load big matrix with "analytical" Coulomb waves
if EnergySpectrum | DoublyDiffSpectrum & ~NumericalEgrid
  load(FileName)
end

% Plot evolution of norm
figure(99)
plot(tVector, normSq, 'k-', 'linewidth', 1.5)
set(gca, 'fontsize', 15)
xlabel('time [a.u.]')
ylabel('|\Psi(t)|^2')
grid on

% Determine and present the angular spectrum - differential in position 
% angle, not asymptotic ejection angle
if AngularSpectrum
  DensMatThetaBefore = 2*dt*h*real(DensMatTheta);  
  EvolveAndAnalyseAngular
end

% Calculate and present singly differential energy spectrum
if EnergySpectrum
  if NumericalEgrid
    EvolveAndAnalyseSingle  
  else
    SingleDiffWithoutInterp
  end
end

% Calculate and present doubly differential spectrum - and singly
% differential spectra in energy and in asymptotic ejection angle
if DoublyDiffSpectrum
  % Assign prefactors
  DensMatDouble = DensMatDouble*dt*sqrt(pi/2)*h^2;
  DoubleDiffWithoutInterp
end
