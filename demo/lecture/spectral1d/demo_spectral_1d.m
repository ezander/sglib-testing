function demo_spectral_1d

clf;
ui = linspace(-1,1)';
plot(ui, G(X(ui))); legend_add('original model');
hold all;

% size of the basis
n = 3;

% Interpolation
u = legendre_rule(n);
V = legendre_eval(n, u);
g = V' \ G(X(u));

plot(ui, g' * legendre_eval(n, ui)); legend_add('interpolation');
plot(u, g' * legendre_eval(n, u), 'rx'); legend_add('interpolation points');

% L2-projection 
[u, w] = legendre_rule(5);
h = legendre_norm(n).^2;
g = legendre_eval(n, u) * diag(w) * G(X(u)) ./ h;

plot(ui, g' * legendre_eval(n, ui)); legend_add('projection');
plot(u, g' * legendre_eval(n, u), 'bx'); legend_add('integration points');

% Galerkin-projection 
L = 2;
xi = [X(0); X(1)-X(0)]; % Poor man's version to get the coefficients of X expanded in the P_i (not good for general code)

triples = legendre_triples(n, L); % triples is a tensor of order 3
Delta = tensor_multiply(triples, xi, 3, 1); % contract tensor over dim 3 and 1 respectively, giving a matrix
rhs = unitvector(1,n);
g = Delta\rhs;

plot(ui, 0.01+g' * legendre_eval(n, ui)); legend_add('Galerkin'); % shift a little bit because otherwise interpolation and Galerking coincide here




function [x, w] = legendre_rule(n)
% LEGENDRE_RULE Return nodes and weights for Gauss-Legendre quadrature
[x, w] = polysys_int_rule('P', n)
w=w';
x=x';

function y=legendre_eval(n, x)
% LEGENDRE_EVAL Evaluate basis of size N of Legendre Polynomials
V=gpcbasis_create('P', 'm', 1, 'p', n-1);
y = gpcbasis_evaluate(V, x');

function h=legendre_norm(n)
% LEGENDRE_NORM Return norm Legendre Polynomials w.r.t. measure induced by U[-1,1]
V=gpcbasis_create('P', 'm', 1, 'p', n-1);
h = gpcbasis_norm(V);

function Delta=legendre_triples(n, l)
% LEGENDRE_TRIPLES Return triples E(P_iP_jP_k) for i,j<n and k<l
V_Y=gpcbasis_create('P', 'm', 1, 'p', n-1);
V_X=gpcbasis_create('P', 'm', 1, 'p', l-1);
Delta = gpcbasis_triples(V_Y, V_Y, V_X);


function y=G(x)
% G The original model 
y = 1./x;

function x=X(u)
% X Transform from a uniform U[-1,1] random variable to the distribution of
% X
a = 0.5;
b = 3;
x = (a+b)/2 + u * (b-a)/2;
