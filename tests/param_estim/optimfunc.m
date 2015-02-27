function [f, df] = optimfunc(x, y, ampl, alpha, omega, phi)
dy = y - ampl * exp(-alpha * x).*sin(omega * x + phi);
f = sum(dy.*dy);
dfdampl = 2 * dot(dy, -exp(-alpha * x).*sin(omega * x + phi));
dfdalpha = 2 * dot(dy, x.*exp(-alpha * x).*sin(omega * x + phi));
dfdomega = 2 * dot(dy, -x.*exp(-alpha * x).*cos(omega * x + phi));
dfdphi = 2 * dot(dy, -exp(-alpha * x).*cos(omega * x + phi));
df = [dfdampl; dfdalpha; dfdomega; dfdphi];
