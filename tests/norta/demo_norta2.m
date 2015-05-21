% Set up the distributions and number samples
N=10000;
X1 = gendist_create('beta', {3, 2});
X2 = gendist_create('lognormal', {1, 0.4});
rho_X = 0.2;

% Generate independent samples
xi_tilde = randn(2, N);

% Generate independend samples and plot
multiplot_init(4, 2, 'ordering', 'row');
x1 = gendist_stdnor(xi_tilde(1,:), X1);
x2 = gendist_stdnor(xi_tilde(2,:), X2);
plot_correlations([x1;x2]);


% Now generate dependend samples with the NORTA method 

% a) Transform the covariance
rho_N = transform_cov_to_normal_fsolve(X1, X2, rho_X, true);
% rho_N = transform_cov_to_normal_fsolve(X1, X2, rho_X, false);
% rho_N = transform_cov_to_normal_hermite(X1, X2, rho_X, 3);
% rho_N = transform_cov_to_normal_hermite(X1, X2, rho_X, 4);
% rho_N = transform_cov_to_normal_hermite(X1, X2, rho_X, 5);

% b) Set up matrix to make the Gaussians have the right covariance
L = gaussian_transform_matrix(rho_N);

% c) Rotate the Gaussians, transform and plot
xi = L * xi_tilde;
x1 = gendist_stdnor(xi(1,:), X1);
x2 = gendist_stdnor(xi(2,:), X2);
plot_correlations([x1;x2]);

% Check the covariance
C_X = cov(x1,x2);
strvarexpand('Covariance should be $rho_X$ and is $C_X(1,2)$.');
