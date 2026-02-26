function F=SpherBess(kr,L)

% The spherical Bessel function j_l(kr).

F=sqrt(pi/2./kr).*besselj(L+1/2,kr);
