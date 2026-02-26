function [DoS LastNeg] = DensityOfStates(EgridRaw)

% For a given energy grid, this function
% calculates the density of states for the
% positive part of the spectrum.

% Number of points
N=length(EgridRaw);

if min(EgridRaw)<0
  LastNeg = max(find(EgridRaw<0));
else
  LastNeg = 0;  
end
EgridRaw =EgridRaw((LastNeg+1):N);
ll=length(EgridRaw);
DoS = zeros(1,ll);
DoS(1) = 1/(EgridRaw(2)-EgridRaw(1));
for nn=2:(ll-1) 
  DoS(nn)=2/(EgridRaw(nn+1)-EgridRaw(nn-1));  
end
DoS(end) =1/(EgridRaw(ll)-EgridRaw(ll-1)); 
DoS=DoS.';