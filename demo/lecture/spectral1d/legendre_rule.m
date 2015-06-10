function [x, w] = legendre_rule(n)
% LEGENDRE_RULE Return nodes and weights for Gauss-Legendre quadrature
[x, w] = polysys_int_rule('P', n);
w=w';
x=x';

