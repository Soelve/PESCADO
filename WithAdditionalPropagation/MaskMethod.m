Aintegral=AfieldInt(E0,w,Ncycle,cep,t);
%PhaseScrinzi=KThMatK.^2/2*t-Aintegral*KThMatK.*cos(KThMatTh);
Phase=KThMatK.^2/2*t + Aintegral*KThMatK.*cos(KThMatTh);
for L=0:lmax
  BesselMat=SpherBess(KRmatK.*KRmatR,L);
  fL=Psi(:,L+1);
  IntegralTerm=trapz(rVector,...
      BesselMat.*rVector.*Gamma.*fL);
  bMatrixMask = bMatrixMask + ...
        i^(-L)*exp(+i*Phase).*YlmMat(KThMatTh,L).*...
        IntegralTerm;      
end
