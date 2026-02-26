%
% Inputs for the grid
%
Rmax = 150;        % Box size
lmax = 12;         % Number of partial waves
N = 750;          % Radial grid points 

%
% Input for the absorber
%
eta = 1e-4;           % Strength
Onset = 100;   
AbsPower = 2;         % Power of the monimial

% Input for the Yukawa potential
Alpha=0.0;

%
% Inputs for the interaction
%
Ncycle = 10;              % # optical cycles
w = .114;                 % Central frequency
E0 = .075;                % Maximum electric field strength
cep = 0;                  % Carrier-envelope phase

% Derivative:
Tpulse = Ncycle*2*pi/w;   % Length of pulse - in time units

% Name of the file with Coulomb waves
FileName = 'PsiMatrixCase400nm.mat';

% Wether angular distribution should be displayed in 3D or 2D
Plot3D = logical(0);

% The theta grid
Ntheta = 200;

%
% Logical inputs
%
% Wether to plott on the fly on not
PlotOrNot =          logical(1);   
% Which spectra to calculate 
AngularSpectrum =    logical(1);
EnergySpectrum =     logical(1);
DoublyDiffSpectrum = logical(1);

% Whether we should use numerical energy grid
NumericalEgrid = logical(0);
% The energy grid
Emin = .01;
Emax = 1;
dE = 2.5e-3;

%
% Numerical inputs for the propagation:
%
% Numerical time step
StepsPerCycle = 1000;
dt = 2*pi/w/StepsPerCycle;

% Dimension of Krylov space for the Arnoldi propagator (fixed):
KrylovDim = 20;           