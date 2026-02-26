% This script calculates the contribution to the energy-spectrum
% for t=T until infinity. Semianalytically

% Cutoff - in imag and real
%CutoffIm = -1e-6;
CutoffIm = -1e-12;
%CutoffIm = 0;
CutoffReal = 0;

% Allocate
disp('Setting up large matrix with couplings for angular distribution')
BigPhiMat = zeros(N,N,lmax+1,lmax+1);
% Set up matrix elements of the CAP between different NH eigenstates
for L = 0:lmax
  for LL = 0:lmax
    %BigPhiMat(:,:,L+1,LL+1) = h*UmatComplex(:,:,L+1)'*...
    %    diag(Gamma)*UmatComplex(:,:,LL+1);
    BigPhiMat(:,:,L+1,LL+1) = h*UmatComplex(:,:,LL+1)'*...
        diag(Gamma)*UmatComplex(:,:,L+1);
  end
end

% Allocate
DensMatThetaAfter = zeros(lmax+1,lmax+1);
% Set up the time-integral of the full contribution beyond t=T
disp('Calculating contribution to angular distribution for t>T.')
for L = 0:lmax
  for LL = 0:lmax
    Entry = 0;
    for n = 1:N
      EnergyN = EnergyComplex(n,L+1);
      %if 1>0
      if imag(EnergyN) < CutoffIm 
      %if imag(EnergyN) < 0 & real(EnergyN) > 0
        for m = 1:N
          EnergyM = EnergyComplex(m,LL+1);
          ProjN = Projections(n,L+1);
          ProjM = Projections(m,LL+1);
          Entry = Entry + ProjN*conj(ProjM)/...
             (EnergyN-conj(EnergyM))*BigPhiMat(m,n,L+1,LL+1);
          % EnergyM = EnergyComplex(m,LL+1);
          % ProjN = Projections(n,L+1);
          % ProjM = Projections(m,LL+1);
          % Entry = Entry + ProjN*conj(ProjM)/...
          %    (EnergyN-conj(EnergyM))*BigPhiMat(m,n,L+1,LL+1);
        end
      end
    end
    DensMatThetaAfter(L+1,LL+1) = Entry;
  end
end

% Prefactor
DensMatThetaAfter = 2*imag(DensMatThetaAfter);

% Add together for total matrix
DensMatThetaTotal = DensMatThetaBefore + DensMatThetaAfter;

% Plot l-distribution of absoroption
figure(77)
%bar(0:lmax, diag(DensMatThetaTotal),'r')
pcolor(0:lmax, 0:lmax, log(abs(DensMatThetaTotal)))
colorbar
xlabel('L')

% Plot angular distribution
% Calculate distribution for matrix
dPdTheta = AngularDist(ThetaGrid,DensMatThetaTotal,lmax);
dPdThetaBefore = AngularDist(ThetaGrid,DensMatThetaBefore,lmax);
dPdThetaAfter = AngularDist(ThetaGrid,DensMatThetaAfter,lmax);
%dPdTheta = dPdThetaBefore + dPdThetaAfter;
FigureNr = 25;
if Plot3D
  %col=[.8 .1 .1];
  col=[0.80,0.53,0.00];
  Handle=PlotAngularDist(dPdTheta,ThetaGrid,FigureNr,col);
else
  figure(FigureNr)
  style = 'k-';  
  %hold on
  polarplot(pi/2-ThetaGrid.', dPdTheta,style)
  hold on
  polarplot(pi/2+ThetaGrid.', dPdTheta,style)
  %polar(pi/2-ThetaGrid.', dPdThetaBefore,'b-.')
  %polar(pi/2+ThetaGrid.', dPdThetaBefore,'b-.')
  %polar(pi/2-ThetaGrid.', dPdThetaAfter,'r--')
  %polar(pi/2+ThetaGrid.', dPdThetaAfter,'r--')
  hold off
end
%legend('Total','Total','Før','Før','Etter', 'Etter')