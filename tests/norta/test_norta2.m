function test_norta

clc
format compact
format short g
%test_normal_decomp
test_cov_decomp



function test_cov_decomp
N=10000;
X1 = gendist_create('beta', {3, 2});
X2 = gendist_create('lognormal', {1, 0.4});
xi1 = gendist_sample(N, X1);
xi2 = gendist_sample(N, X2);
multiplot_init(4);
multiplot; plot_density(X1); hold all; plot_density(xi1, 'type', 'kde');
multiplot; plot_density(X2); hold all; plot_density(xi2, 'type', 'kde');
%multiplot_adjust_range

%rho_X = 0.3;
%rho_N = transform_cov_to_normal(X1, X2, rho_X);
rho_N = 0.4;
rho_X = transform_cov_from_normal(X1, X2, rho_N);


function rho_N = transform_cov_to_normal(X1, X2, rho_X);
assert(false);

function rho_X = transform_cov_from_normal(X1, X2, rho_N);
C = [1, rho_N; rho_N, 1];
L = chol(C, 'lower');
xi = randn(2, 100000);
xi = normal_invcdf(halton_sequence(100000, 2)');

[x1, w1] = gauss_hermite_rule(7);
[x, w] = tensor_mesh({x1, x1}, {w1, w1});

xc = L*x;
assert( norm(C - xc * diag(w) * xc') < 1e-10 )

mu1 = gendist_moments(X1);
mu2 = gendist_moments(X2);
xx * diag(w) * xx'




%assert( norm(C - cov((L*xi)')) < 0.01 )

xi2 = L * xi;
eta1 = [gendist_stdnor(xi(1,:), X1); gendist_stdnor(xi(2,:), X2)];
eta2 = [gendist_stdnor(xi2(1,:), X1); gendist_stdnor(xi2(2,:), X2)];

multiplot; 
plot_density(eta1(1,:)); hold all; plot_density(eta2(1,:))
multiplot; 
plot_density(eta1(2,:)); hold all; plot_density(eta2(2,:))
multiplot_adjust_range('separate', 'rows')

%%
clc

mu1 = gendist_moments(X1);
mu2 = gendist_moments(X2);
xx = [gendist_stdnor(x(1,:), X1)-mu1; gendist_stdnor(x(2,:), X2)-mu2];
cov(eta1')
xx * diag(w) * xx'
xx = [gendist_stdnor(xc(1,:), X1)-mu1; gendist_stdnor(xc(2,:), X2)-mu2];
xx * diag(w) * xx'
cov(eta2')


%x1 = gendist_stdnor(



function test_normal_decomp
m = 3;

% Set up a covariance matrix
A = rand(m, m);
C = A*A';
C = 0.5 * (C + C');
assert(norm(C-C')<1e-10) % make sure C is symmetric
assert(all(eig(C))>0) % make sure C is positive

% Do a Cholesky factorisation
L = chol(C, 'lower'); % need to use 'lower' (default is right upper matrix)
assert( norm(C - L*L' )<1e-10); % make sure L really is a factorisation of C

%
N = 10000;
xi = randn(m, N);
g = L * xi;

C_g = cov(g');
C
C_g
rel_err = norm(C_g-C)/norm(C)




