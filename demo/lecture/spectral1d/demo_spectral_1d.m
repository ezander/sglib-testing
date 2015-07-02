function demo_spectral_1d

clf;
ui = linspace(-1,1)';
plot(ui, G(X(ui))); legend_add('original model');
hold all;

% size of the basis
n = 4;

% Interpolation
u = legendre_rule(n); % interpolation points
% u = 2*rand(n,1)-1 % this would also do
V = legendre_eval(n, u); % compute the Vandermonde matrix
g = V' \ G(X(u)); 
g_I = g;

plot(ui, g' * legendre_eval(n, ui)); legend_add('interpolation');
plot(u, g' * legendre_eval(n, u), 'rx'); legend_add('interpolation points');

% L2-projection 
[u, w] = legendre_rule(5);
h = legendre_norm(n).^2;
g = legendre_eval(n, u) * diag(w) * G(X(u)) ./ h;
g_P = g;

plot(ui, g' * legendre_eval(n, ui)); legend_add('projection');
plot(u, g' * legendre_eval(n, u), 'bx'); legend_add('integration points');

% Galerkin-projection 
L = 4;
[u, w] = legendre_rule(5);
h = legendre_norm(L).^2;
xi = legendre_eval(L, u) * diag(w) * X(u) ./ h;

triples = legendre_triples(n, L); % triples is a tensor of order 3
Delta = tensor_multiply(triples, xi, 3, 1); % contract tensor over dim 3 and 1 respectively, giving a matrix
rhs = unitvector(1,n);
g = Delta\rhs;
g_G = g;

plot(ui, g' * legendre_eval(n, ui)); legend_add('Galerkin'); % shift a little bit because otherwise interpolation and Galerking coincide here

norm(g_G - g_I)

hold off

function y=G(x)
% G The original model 
y = 1./x;

function x=X(u)
% X Transform from a uniform U[-1,1] random variable to the distribution of
% X
a = 0.5;
b = 3;
x = (a+b)/2 + u * (b-a)/2;
x = (a+b)/2 + u * (b-a)/2 + 0.3 * u.^2;
