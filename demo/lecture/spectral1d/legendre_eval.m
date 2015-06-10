function y=legendre_eval(n, x)
% LEGENDRE_EVAL Evaluate basis of size N of Legendre Polynomials
V=gpcbasis_create('P', 'm', 1, 'p', n-1);
y = gpcbasis_evaluate(V, x');

