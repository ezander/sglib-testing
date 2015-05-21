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


