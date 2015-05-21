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

rho_N = 0.4;

rho_X = transform_cov_from_normal_quad(X1, X2, rho_N)
rho_N1 = transform_cov_to_normal_fsolve(X1, X2, rho_X, true)

rho_X = transform_cov_from_normal_qmc(X1, X2, rho_N)
rho_N1 = transform_cov_to_normal_fsolve(X1, X2, rho_X, false)

%p = setup_polynomial(X1, X2, 10, rho_X)
rho_N1 = transform_cov_to_normal_hermite(X1, X2, rho_X, 3)
rho_N1 = transform_cov_to_normal_hermite(X1, X2, rho_X, 4)
rho_N1 = transform_cov_to_normal_hermite(X1, X2, rho_X, 5)


function rho_N = transform_cov_to_normal_fsolve(X1, X2, rho_X, quad)
if quad
    func = @(rho_N)(transform_cov_from_normal_quad(X1, X2, rho_N)-rho_X);
else
    func = @(rho_N)(transform_cov_from_normal_qmc(X1, X2, rho_N)-rho_X);
end
opts = optimset('Display', 'off');
[rho_N, ~, flag] = fsolve(func, rho_X, opts);
if flag~=1
    warning('fsolve had some problems');
end

function rho_N = transform_cov_to_normal_hermite(X1, X2, rho_X, M)
p = setup_polynomial(X1, X2, M, rho_X);
rs = roots(p);
ind=(rs>=-1.1) & (rs<=1.1) & (imag(rs)==0);
rho_N=unique(rs(ind));
if isempty(rho_N)
    error( 'the inverse of the transform polynomial could not be found.' );
elseif length(rho_N)>1
    error( 'the inverse of the transform polynomial is not unique.' );
end


function p_i = setup_polynomial(X1, X2, M, rho_X)
%[xi, w] = gauss_hermite_rule(min(M,7));

deg = min(M, 7);
% Compute the coefficients a_i
a_i = hermite_expand(X1, M, deg);
b_i = hermite_expand(X2, M, deg);
p_i = a_i .* b_i .* factorial((0:M)');
p_i(1) = p_i(1) - rho_X;
p_i(1) = - rho_X;
p_i = p_i(end:-1:1)';

function a_i = hermite_expand(X, M, deg)
V=gpcbasis_create('H', 'm', 1, 'p', M);
[xi,w]=gpc_integrate([], V, deg);
a_i = binfun(@times, gpcbasis_evaluate(V, xi), gendist_stdnor(xi, X)) * w ./ gpcbasis_norm(V, 'sqrt', false);

function rho_X = transform_cov_from_normal_quad(X1, X2, rho_N)
% TRANSFORM_COV_FROM_NORMAL_QUAD Compute covariance from Gaussian cov using quadrature.

% Compute Choleski decomposition of covariance matrix
C = [1, rho_N; rho_N, 1];
L = chol(C, 'lower');

% Generate Gauss-Hermite integration points
[xi1, w1] = gauss_hermite_rule(7);
[xi, w] = tensor_mesh({xi1, xi1}, {w1, w1});

% Compute the points for correlated Gaussians
% (and make sure that they have the correct correlation structure)
g = L*xi;
assert( norm(C - g * diag(w) * g') < 1e-10 )

% Compute mean and transformed points for distribution X1
mu1 = gendist_moments(X1);
x1 = gendist_stdnor(g(1,:), X1); 

% Compute mean and transformed points for distribution X2
mu2 = gendist_moments(X2);
x2 = gendist_stdnor(g(2,:), X2);

% Compute covariance cov(X1, X2) using the integration rule 
rho_X = ((x1-mu1).*(x2-mu2))*w;


function rho_X = transform_cov_from_normal_qmc(X1, X2, rho_N)
% TRANSFORM_COV_FROM_NORMAL_QMC Compute covariance from Gaussian cov using quasi Monte Carlo.

% Compute Choleski decomposition of covariance matrix
C = [1, rho_N; rho_N, 1];
L = chol(C, 'lower');

% Generate QMC points
N = 100000;
xi = normal_invcdf(halton_sequence(N, 2)');

% Compute the points for correlated Gaussians
% (and make sure that they have approx. the correct correlation structure)
g = L*xi;
assert( norm(C - cov(g')) < 1e-3 )

% Compute mean and transformed points for distribution X1
mu1 = gendist_moments(X1);
x1 = gendist_stdnor(g(1,:), X1); 

% Compute mean and transformed points for distribution X2
mu2 = gendist_moments(X2);
x2 = gendist_stdnor(g(2,:), X2);

% Compute covariance cov(X1, X2) using QMC
rho_X = (x1-mu1)*(x2-mu2)'/N;
