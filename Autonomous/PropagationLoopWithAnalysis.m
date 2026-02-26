% Initiate time
t=0;
nrr=1;

% Plot wave function
if PlotOrNot    
  figure(2)
  %
  % Plot field and progress
  %
  subplot(5,1,1)
  ttt=0:dt:Tpulse; 
  plot(ttt,Afield(E0,w,Ncycle,cep,ttt))
  hold on
  plot(ttt,Efield(E0,w,Ncycle,cep,ttt))
  PlotPoint=plot(t,0,'r*');
  hold off
  xlabel('t')
  subplot(5,1,2)
  %
  % Plot |f_{l,m}(r,t)|^2 for l=0, 1 og 2.
  %
  l0plott=plot(rVector,abs(Psi(:,1)).^2,'b-');
  title('l=0')
  xlabel('r')
  subplot(5,1,3)
  l1plott=plot(rVector,abs(Psi(:,2)).^2,'r-');
  title('l=1')
  xlabel('r')
  subplot(5,1,4)
  l2plott=plot(rVector,abs(Psi(:,3)).^2,'g-');
  title('l=2')
  xlabel('r')
  subplot(5,1,5)
  l3plott=plot(rVector,abs(Psi(:,4)).^2,'m-');
  title('l=3')
  xlabel('r')
end

% Construct theta-grid
ThetaGrid=linspace(0,pi,Ntheta);
dPdTheta=zeros(1,Ntheta);
SpherHarmMat=zeros(Ntheta,lmax+1);

% Construct the Spherical Harmonics
for ll=0:lmax;
  LegVector=legendre(ll,cos(ThetaGrid)); % Legendre polynomial
  LegVector=LegVector(1,:).';
  % Include normalization
  SpherHarmMat(:,ll+1)=sqrt((2*ll+1)/(4*pi))*LegVector; 
end

%
% Initiate and allocate density matrices for spectra
%
if AngularSpectrum
  DensMatTheta=zeros(lmax+1,lmax+1);
end

if EnergySpectrum
  % DensMatEnergy=zeros(N,N,lmax+1);
  DensMatEnergy = zeros(N-Nt+1,N, lmax+1);
end

if DoublyDiffSpectrum
  DensMatDouble=zeros(N,N,lmax+1,lmax+1);
end

% 
% End of allocation
%

% Propagation loop - with field
counter = 0;
while t < Tpulse 
  % Half step with absorber  
  Psi=diag(exp(-dt/2*Gamma))*Psi;
  % Propagate Hermitian part, Arnoldi version
  Psi=ArnoldiProp(H0rad,Sr,Sl,rVector,Pz,...
  AngularOp,AngularOp2,E0,w,Ncycle,cep,t+dt/2,...
  Psi,KrylovDim,dt,h);
  % Another half step with absorber
  Psi=diag(exp(-dt/2*Gamma))*Psi;
  
  t=t+dt;                           % Update time
  counter = counter+1;
  % Update plots - progress and f_l for l=0 and 1 (if requested)
  if PlotOrNot & mod(counter,100)==0
    set(PlotPoint,'xdata',t)
    drawnow
    set(l0plott,'ydata',abs(Psi(:,1)).^2)
    set(l1plott,'ydata',abs(Psi(:,2)).^2)
    set(l2plott,'ydata',abs(Psi(:,3)).^2)
    set(l3plott,'ydata',abs(Psi(:,4)).^2)
    drawnow
    %pause
  elseif ~PlotOrNot & mod(counter,100)==0
    % Write progress to screen    
    disp(['Progress: ', num2str(t/Tpulse*100),' %'])
  end
  
  % Calculate and store norm
  normSq(nrr)=NormMatrixPsi(Psi,h)^2;       
  tVector(nrr)=t;                           % Vector with time grid
  nrr=nrr+1;
  
  %
  % Calculate density matrices for spectra
  %
  
  % Angular spectrum only - by direct measurement
  if AngularSpectrum
    DensMatTheta = DensMatTheta + Psi'*diag(Gamma)*Psi;   
  end
  
  % Energy spectrum only
  if EnergySpectrum
    for L=0:lmax
      DensMatEnergy(:,:,L+1)=DensMatEnergy(:,:,L+1)+...
          Psi(Nt:N,L+1)*Psi(:,L+1)';
    end
  end
  
  % Fully differential distribution - by the PESCADO method
  if DoublyDiffSpectrum
    for L=0:lmax
      for Lp=0:lmax
        DensMatDouble(:,:,L+1,Lp+1) = ...
        DensMatDouble(:,:,L+1,Lp+1)+Psi(:,L+1)*Psi(:,Lp+1)';  
      end  
    end
  end       % End if-statement  
end         % End of time-loop


%
% Half step back - to clean up
%
% Angular spectrum only - by direct measurement
if AngularSpectrum
  DensMatTheta = DensMatTheta - 0.5*Psi'*diag(Gamma)*Psi;   
end
% Energy spectrum only
if EnergySpectrum
  for L=0:lmax
    DensMatEnergy(:,:,L+1)=DensMatEnergy(:,:,L+1) - ...
        0.5*Psi(Nt:N,L+1)*Psi(:,L+1)';
  end
end
% Fully differential distribution - by the PESCADO method
if DoublyDiffSpectrum
  for L=0:lmax
    for Lp=0:lmax
      DensMatDouble(:,:,L+1,Lp+1) = ...
      DensMatDouble(:,:,L+1,Lp+1) - 0.5*Psi(:,L+1)*Psi(:,Lp+1)';  
    end  
  end
end       % End if-statement  

% Fill with zeros - for bookeeping
if EnergySpectrum
  aux = zeros(N,N);
  NewDensMatEnergy = zeros(N,N,lmax+1);
  for L=0:lmax
    NewDensMatEnergy(:,:,L+1) = aux;
    NewDensMatEnergy(Nt:N,:,L+1) = DensMatEnergy(:,:,L+1);
  end
  DensMatEnergy = NewDensMatEnergy;
 clear NewDensMatEnergy
end

return
  
% Correct with dt - and copy eta(t=T)
%if DoublyDiffSpectrum
%  DensMatDouble = DensMatDouble*dt; 
%  DensMatDoubleBefore = DensMatDouble;
%end

% Find spectrum of Hermitian Hamiltonian
CalculateSpectra

% Find spectrum of Non-Hermitian Hamiltonian
CalculateNonHermSpectra

% Ionization probability differential in energy (not interpolated)
if EnergySpectrum
  EdistRaw=zeros(N,lmax+1);
  for L=0:lmax
    for EnInd=1:N
      EdistRaw(EnInd,L+1)=real(Umat(:,EnInd,L+1)'*diag(Gamma)*...
      DensMatEnergy(:,:,L+1)*Umat(:,EnInd,L+1));
    end
  end
  % Proper prefactors
  EdistRaw=EdistRaw*2*h^2*dt;
  % Interpolation
  InterpolateEdist
end

%PsiT = Psi;          % Copy state

% Calculate energy spectra obtained for t <= T
%PresentResults

if AngularSpectrum
 DensMatThetaBefore = 2*dt*h*real(DensMatTheta);  % Introduce proper factors
 EvolveAndAnalyseAngular
end

% Calculate additional energy contribution from after the pulse
if EnergySpectrum
  EvolveAndAnalyseSingle
end