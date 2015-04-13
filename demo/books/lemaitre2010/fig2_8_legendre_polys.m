function fig2_8_legendre_polys
% Fig 2.8 page 36 (lemaitre & knio spectral methods)

V = gpcbasis_create('P', 'm', 1, 'p', 6);

x = linspace(-1, 1);
y = gpcbasis_evaluate(V, x);

clf
plot(x, y);
set(gcf, 'defaulttextinterpreter', 'latex');

xlabel('$\xi$');
ylabel('$P_i(\xi)$');
title('Legendre polynomials');
