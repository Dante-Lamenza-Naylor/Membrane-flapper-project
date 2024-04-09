function sol = Theodorsen(s)
% s is reduced frequency pi * c * f / U_inf
sol = besselk(1,s*1i)/(besselk(1,s*1i) + besselk(0,s*1i));