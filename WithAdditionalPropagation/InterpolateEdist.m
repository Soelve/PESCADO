% This script calculates the energy distribution 
% which emerges from the absorption on an 
% energy grid which is the same for all partial waves

% Set up grid for interpolation
dE=.0025;
Evector=Emin:dE:Emax;
EdistInterpolated=zeros(length(Evector),lmax+1);

for L=0:lmax
  EgridRaw=Emat(:,L+1);  
%  if min(EgridRaw)<0
%    LastNeg = max(find(EgridRaw<0));
%  else
%    LastNeg = 0;  
%  end
%  EgridRaw =EgridRaw((LastNeg+1):N);
  [DoS LastNeg] = DensityOfStates(EgridRaw);
  EdistRawL=EdistRaw((LastNeg+1):N,L+1);
%  ll=length(EgridRaw);
%  DoS = zeros(1,ll);
%  DoS(1) = 1/(EgridRaw(2)-EgridRaw(1));
%  for nn=2:(ll-1) 
%    DoS(nn)=2/(EgridRaw(nn+1)-EgridRaw(nn-1));  
%  end
%  DoS(end) =1/(EgridRaw(ll)-EgridRaw(ll-1)); DoS=DoS.';
  EdistInterpolated(:,L+1)=spline(EgridRaw((LastNeg+1):N),...
      EdistRawL.*DoS,Evector);
end