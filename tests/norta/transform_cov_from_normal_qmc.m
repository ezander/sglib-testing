function rho_X = transform_cov_from_normal_qmc(X1, X2, rho_N)
% TRANSFORM_COV_FROM_NORMAL_QMC Compute covariance from Gaussian cov using quasi Monte Carlo.

% Compute Choleski decomposition of covariance matrix
C = [1, rho_N; rho_N, 1];
[L, n] = covariance_decomp(C);

% Generate QMC points
N = 100000;
xi = normal_samples(n, N, true);

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
