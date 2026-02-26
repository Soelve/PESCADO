% Create 3d plot of the angular distrubution
% Grid in r, theta and phi
R=dPdTheta;
Th=ThetaGrid;
Phi=0:.05:(2*pi+.05);

% Sizes of the grids
Nth=length(Th);
Nphi=length(Phi);

% Convert to Cartesian coordinates
for ll=1:Nth
  rr=R(ll);
  for kk=1:Nphi
    X(ll,kk)=rr*sin(Th(ll))*cos(Phi(kk));
    Y(ll,kk)=rr*sin(Th(ll))*sin(Phi(kk));
    Z(ll,kk)=rr*cos(Th(ll));
  end
end

% Make surface plot
figure(7)
PlotHandle=surf(X,Y,Z);
shading interp
set(PlotHandle,'FaceColor',[.3 .5 .2],'FaceAlpha',.8)   % Set colour and transparency
% Adjust ligthing
lighting flat
camlight('headlinght')
camlight('left')
% Set axes equal and remove them
axis off
axis equal
