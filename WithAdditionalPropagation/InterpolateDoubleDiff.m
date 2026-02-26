% This script calculates the energy distribution 
% which emerges from the absorption on an 
% energy grid which is the same for all partial waves

% Energy vectors
Evector=Emin:dE:Emax;
Nenergy=length(Evector);

% Allocate
DoublyDiffDist = zeros(Nenergy,Ntheta);
Tot = zeros(1,Nenergy);

% Grid
[EEx EEy]=meshgrid(Evector,Evector);

for L=0:lmax
  sigmaL = angle(gammaZ(L+1-1i./sqrt(2*Evector)));
  EgridRawL=Emat(:,L+1);  
  [DoSL LastNegL] = DensityOfStates(EgridRawL); 
  for Lp=0:lmax
    disp([L Lp])
    clear Aux1 Aux2 Aux3
    sigmaLp = angle(gammaZ(Lp+1-1i./sqrt(2*Evector)));
    EgridRawLp=Emat(:,Lp+1);  
    [DoSLp LastNegLp] = DensityOfStates(EgridRawLp); 
    [RawL RawLp]=meshgrid(EgridRawL((LastNegL+1):N),...
        EgridRawLp((LastNegLp+1):N));

    Aux1=h^2*Umat(:,(LastNegL+1):N,L+1)'*diag(Gamma)*...
        DensMatDouble(:,:,L+1,Lp+1)*Umat(:,(LastNegLp+1):N,Lp+1);
    
    Aux2 = diag(sqrt(DoSL))*Aux1*diag(sqrt(DoSLp));
    Aux3 = interp2(RawL,RawLp,Aux2.',EEx,EEy,'spline');    

    ProjectionPart=diag(Aux3);
    Tot=Tot+ProjectionPart.';
    % Product of spherical harmonics
    AngleFunc=Ylm(ThetaGrid,L).*conj(Ylm(ThetaGrid,Lp));
    % Augment d^2P/dE d\theta with the term of l and l'
    DoublyDiffDist = DoublyDiffDist + ...
        2*real(i^(Lp-L)*exp(1i*(sigmaL-sigmaLp)).'.*...
        ProjectionPart*AngleFunc);
  end
end

figure(23)
[EmatMesh ThetaMat]=meshgrid(Evector,ThetaGrid);
pcolor(EmatMesh.*sin(ThetaMat),EmatMesh.*cos(ThetaMat),...
    (DoublyDiffDist.'>0).*DoublyDiffDist.')
shading interp
axis([0 .5 -1 1])
daspect([1 1 1])
C=colormap; C(1,:)=[1 1 1]; colormap(C)
set(gca,'fontsize',15)
xlabel('E cos \theta (a.u.)')
ylabel('E sin \theta (a.u.)')