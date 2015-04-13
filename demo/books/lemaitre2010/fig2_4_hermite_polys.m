function fig2_4_hermite_polys
% Fig 2.4 page 32 (lemaitre & knio spectral methods)

V = gpcbasis_create('H', 'm', 1, 'p', 6);

x = linspace(-3, 3);
y = gpcbasis_evaluate(V, x);

clf
plot(x, y);
set(gcf, 'defaulttextinterpreter', 'latex');

ylim([-20, 25]);
xlabel('$\xi$');
ylabel('$H_i(\xi)$');
title('Hermite polynomials');
