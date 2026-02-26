% Plot the norm of the wave function as a function of time
figure(3)
hold off
plot(tVektor,normSq,'--', 'color', [1 0 0], 'linewidth',1.)
title(['E_0= ',num2str(E0),', \omega= ',num2str(w),', \Delta t=',num2str(dt)])
xline(Tpulse,'k--');
%xline(2*Tpulse+5,'k--');
xlabel('Time [a.u.]')
ylabel('|\Psi(t)|^2')
set(gca,'fontsize',15)
grid on

% Calculate dP/d \theta
if AngularSpectrum
  dPdTheta=AngularDist(ThetaGrid,DensMatTheta,lmax);
  dPdTheta1=dPdTheta;
  FigNr=111;
  col=[.3 .5 .2];
  Handle=PlotAngularDist(dPdTheta,ThetaGrid,FigNr,col);
  
  % Plot l-distribution (including correlation elements)
  figure(12)
  pcolor(0:lmax,0:lmax,abs(DensMatTheta))
end

%
% Plot Energy-distribution
%

if EnergySpectrum
  figure(13)
  hold off
  % Plot total dP/dE (sum over L)
  EdistTotal=sum(EdistInterpolated,2);
  %semilogy(Evector,EdistTotal,'g-.','linewidth',1.2)
  hold on
  plot(Evector,EdistTotal,'g-.','linewidth',1.2)
  % Calculate negative probability - and write to screen
  NegativeContributionToEnergySpectrum=- ...
      trapz(Evector,(EdistTotal<0).*EdistTotal)

  Vax=axis;
  axis([0 1.5 -.2 Vax(4)])
  xlabel('Energy (a.u.)')
  ylabel('dP/dE (a.u.)')
  set(gca,'fontsize',15)
  grid on
end

if DoublyDiffSpectrum
  InterpolateDoubleDiff
  %
  % Plot dP/dE and dP/dOmega from d^2P / dE dOmega
  %
  Edist2=trapz(ThetaGrid,sin(ThetaGrid).*DoublyDiffDist,2)*2*pi;
  figure(15)
  hold off
  semilogy(Evector,Edist2,'--')
  hold off
  title('PESCADO')
  xlabel('Energy (a.u.)')
  ylabel('dP/dE')
  grid on
  
  % Integrate out energy dependence
  dPdTheta2=trapz(Evector,DoublyDiffDist,1);
  
  FigNr=18;
  col=[.2 .7 .1];
  Handle=PlotAngularDist(dPdTheta2,ThetaGrid,FigNr,col);
  title('PESCADO')
end

if MaskMeth
  % Energy distribtuion
  EdistMask = trapz(ThetaGrid,...
      2*pi*sin(ThetaGrid).'.*abs(bMatrixMask).^2.*kVector);
  figure(25)
  hold off
  semilogy(kVector.^2/2,EdistMask,'r--')
  title('Mask method')
  xlabel('Energy (a.u.)')
  ylabel('dP/dE')
  grid on

  % Angular distribution
  dPdThetaMask=trapz(kVector,abs(bMatrixMask).^2*diag(kVector.^2),2);
  %dPdThetaMask=trapz(kVector,abs(bMatrixMask).^2,2);
  
  FigNr=27;
  col=[.7 .2 .2];
  Handle=PlotAngularDist(dPdThetaMask,ThetaGrid,FigNr,col);
  title('Mask method')
end

grid on
set(gca, 'fontsize', 15)