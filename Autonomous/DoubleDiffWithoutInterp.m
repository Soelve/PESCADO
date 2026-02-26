% This script calculates the ionization probability singly differential
% in energy. It does so using a precalculated set of Coulomb waves. These
% waves, in turn, are obtained using the Special Functions in Physics
% MATLAB Toolbox.

% Load array with analytical Coulomb waves
%load PsiMatrixCase1.mat
ImagCutoff = 1e-12;

% Allocate
Nenergy_an = length(Egrid);
DoublyDiffDist_Before = zeros(Nenergy_an,Ntheta);
DoublyDiffDist_After = zeros(Nenergy_an,Ntheta);
ProjectionPart = zeros(Nenergy_an,1);

% Assign prefactors
%DensMatDouble = DensMatDouble*dt*sqrt(pi/2)*h^2;

% 2D energy grid
[EEx EEy]=meshgrid(Egrid,Egrid);

% The part obtained for t<=T
disp('Calculating contribution to the doubly differential during int.')
for L=0:lmax
  sigmaL = angle(gammaZ(L+1-1i./sqrt(2*Egrid)));
  PsiMatForL = BigPsiMat(:,:,L+1);
  for Lp=0:lmax
    clear Aux1;
    sigmaLp = angle(gammaZ(Lp+1-1i./sqrt(2*Egrid)));
    PsiMatForLp = BigPsiMat(:,:,Lp+1);
    % Calculate "sandwitch"
    Aux1=PsiMatForL'*diag(Gamma)*...
        DensMatDouble(:,:,L+1,Lp+1)*PsiMatForLp;
    % Extract diagonal and set correct prefactors
    ProjectionPart = diag(Aux1);
    ProjectionPart = ProjectionPart./sqrt(2*Egrid.');
    % Product of spherical harmonics
    AngleFunc=Ylm(ThetaGrid,L).*conj(Ylm(ThetaGrid,Lp));
    % Augment d^2P/dE d\theta with the term of l and l'
    % NB: Necessary to remove factor 2. I dunno why
    DoublyDiffDist_Before = DoublyDiffDist_Before + ...
        real(i^(Lp-L)*exp(1i*(sigmaL-sigmaLp)).'.*...
        ProjectionPart*AngleFunc);
  end  % Lp
end  % L

%
% Contribution for t>T
%
disp('Constructing energydistribution for t>T.')

% Matrix with, eh, matrix elements
disp('Setting up coupling matrices.')
for L = 0:lmax;
  BigMatWithGamma_An(:,:,L+1) = BigPsiMat(:,:,L+1)'*diag(Gamma)*...
      UmatComplex(:,:,L+1);
  BigMatNoGamma_An(:,:,L+1) = BigPsiMat(:,:,L+1)'*UmatComplex(:,:,L+1);
end

% Put it all together into an energy distribution
tic
for L = 0:lmax
  sigmaL = angle(gammaZ(L+1-1i./sqrt(2*Egrid)));
  PsiMatForL = BigPsiMat(:,:,L+1);
  for Lp = 0:lmax    
    disp(['Calculating contribution for L=',num2str(L),...
        ', Lp=',num2str(Lp)'.'])
    sigmaLp = angle(gammaZ(Lp+1-1i./sqrt(2*Egrid)));
    PsiMatForLp = BigPsiMat(:,:,Lp+1);
    % Reset variable
    ProjectionPart = zeros(Nenergy_an, 1);
    % Loop over energy grid (preset)
    for EnergyIndX = 1:Nenergy_an
      Term = 0;         % Initiate relevant element
      % Loop over complex energy N
      for EnergyN = 1:N
        % Check that it is admissible
        if imag(EnergyComplex(EnergyN, L+1)) < -ImagCutoff
          % Loop over complex energy M
          for EnergyM = 1:N
            Term = Term + ...
            Projections(EnergyN,L+1)*conj(Projections(EnergyM,Lp+1))/...
            (EnergyComplex(EnergyN,L+1)-...
            conj(EnergyComplex(EnergyM,Lp+1))) * ...
            BigMatWithGamma_An(EnergyIndX,EnergyN,L+1) * ...
            conj(BigMatNoGamma_An(EnergyIndX,EnergyM,Lp+1));
          end % EnergyM
        end % for if-statement
      end % EnergyN
      % Distribution
      ProjectionPart(EnergyIndX) = Term;
    end % EnergyIndX
    %ProjectionPart = sqrt(pi/2)*h^2*imag(ProjectionPart)./sqrt(2*Egrid.');
    ProjectionPart = sqrt(pi/2)*h^2*ProjectionPart./sqrt(2*Egrid.');
    % Product of spherical harmonics
    AngleFunc=Ylm(ThetaGrid,L).*conj(Ylm(ThetaGrid,Lp));
    % Augment d^2P/dE d\theta with the term of l and l'
    % NB: Necessary to remove factor 2. I dunno why
    DoublyDiffDist_After = DoublyDiffDist_After + ...
        imag(i^(Lp-L)*exp(1i*(sigmaL-sigmaLp)).'.*...
        ProjectionPart*AngleFunc);
  end % Lp
  % Interpolated distribution
end % L
toc

% Add contributions
DoublyDiffDist = DoublyDiffDist_Before + DoublyDiffDist_After;

% Plot doubly differential distribution
figure(27)
[EmatMesh ThetaMat]=meshgrid(Egrid,ThetaGrid);
pcolor(EmatMesh.*sin(ThetaMat),EmatMesh.*cos(ThetaMat),...
    (DoublyDiffDist.'>0).*DoublyDiffDist.')
shading interp
axis([0 .5 -1 1])
daspect([1 1 1])
colormap('hot')
C=colormap; C = flipud(C);
C(1,:)=[1 1 1]; colormap(C)
set(gca,'fontsize',15)
xlabel('E cos \theta (a.u.)')
ylabel('E sin \theta (a.u.)')

%
% Plot energy distribution
%
% Integrate out ejection angle
Edist2=trapz(ThetaGrid,sin(ThetaGrid).*DoublyDiffDist,2)*2*pi;
figure(22)
plot(Egrid, Edist2, 'k-', 'linewidth',1.7)
grid on
xlabel('E [a.u.]')
ylabel('dP_{ion}/dE')
set(gca, 'fontsize', 15)

%
% Plot angular distribution
%
% Integrate out energy
dPdThetaK = trapz(Egrid,DoublyDiffDist,1);
FigureNr = 18;
if Plot3D
  col=[.2 .7 .1];
  Handle=PlotAngularDist(dPdThetaK,ThetaGrid,FigureNr,col);
else
  figure(FigureNr)
  polarplot(pi/2-ThetaGrid,dPdThetaK,'k-')
  hold on
  polarplot(pi/2+ThetaGrid,dPdThetaK,'k-')
  hold off
end