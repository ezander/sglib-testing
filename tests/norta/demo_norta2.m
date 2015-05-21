% Set up the distributions and number samples
N=20000;
X1 = gendist_create('beta', {3, 2});
X2 = gendist_create('lognormal', {1, 0.4});
rho_X = 0.2;

% Generate independent samples
xi_tilde = normal_samples(2, N, true);

% Generate independend samples and plot
multiplot_init(4, 4, 'ordering', 'row');
x1 = gendist_stdnor(xi_tilde(1,:), X1);
x2 = gendist_stdnor(xi_tilde(2,:), X2);
plot_correlations(xi_tilde, 1, 1);
plot_correlations([x1;x2], 1, 3);


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
plot_correlations(xi, 3, 1);
plot_correlations([x1;x2], 3, 3);

% Check the covariance
C_X = cov(x1,x2);
strvarexpand('Covariance should be $rho_X$ and is $C_X(1,2)$.');

