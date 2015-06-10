function h=legendre_norm(n)
% LEGENDRE_NORM Return norm Legendre Polynomials w.r.t. measure induced by U[-1,1]
V=gpcbasis_create('P', 'm', 1, 'p', n-1);
h = gpcbasis_norm(V);

