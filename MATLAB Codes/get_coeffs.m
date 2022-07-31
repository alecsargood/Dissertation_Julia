function [ak, bk, gk, dk, chik] = get_coeffs(k, Du, Dv, u, v, a, b, LI, GM1, GM2)

% λ^2 + ak λ + bκ + (gk λ + dk)E + chik E^2

if LI == 1

ak = ((Du+Dv)*k^2*pi^2) + (u^2) + (4*u*v) +1;

bk = ((Dv*pi^2*k^2) + u^2)*((Du*(pi^2)*(k^2)) + (4*u*v) + 1) - (4*u^3*v);

gk = -6*u*v;

dk = -6*Dv*u*v*(k^2)*(pi^2);

chik = 0;

elseif GM1 == 1

ak = ((Du+Dv)*k^2*pi^2) + b + 1;

bk = ((Du*pi^2*k^2) + b)*((Dv*(pi^2)*(k^2)) + 1);

gk = -2*u/v;

dk = -2*u/v*(Dv*(k^2)*(pi^2) + 1);

chik = 2*u^3 / v^2;

elseif GM2 == 1

ak = ((Du+Dv)*k^2*pi^2) + b + 1;

bk = ((Du*pi^2*k^2) + b)*((Dv*(pi^2)*(k^2)) + 1);

gk = -2*u/v;

dk = -2*u/v*(Dv*(k^2)*(pi^2) + 1) + 2*u^3 / v^2;

chik = 0;

end






end