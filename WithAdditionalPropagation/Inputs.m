%
% Inputs for the grid
%
Rmax = 150;       % Box size
lmax = 12;        % Number of partial waves
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
E0 = .075;                  % Maximum electric field strength
w = .114;                   % Central frequency
cep = 0;                  % Carrier-envelope phase
Textra = 10;             % Propagation after pluse 

% Derivative:
Tpulse = Ncycle*2*pi/w;   % Length of pulse - in time units

% The energy grid
Emin = .001;
Emax = 1;
dE = 1e-2;
dK = .0025;

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
DoublyDiffSpectrum = logical(0);
MaskMeth =           logical(1);

%
% Numerical inputs for the propagation:
%
% Numerical time step
StepsPerCycle = 1000;
dt = 2*pi/w/StepsPerCycle;
% Dimension of Krylov space for the Arnoldi propagator(fixed):
KrylovDim=20;           