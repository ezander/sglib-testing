function rho_X = transform_cov_from_normal_quad(X1, X2, rho_N)
% TRANSFORM_COV_FROM_NORMAL_QUAD Compute covariance from Gaussian cov using quadrature.

% Compute Choleski decomposition of covariance matrix
C = [1, rho_N; rho_N, 1];
L = covariance_decomp(C, 'fillup', true);

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


