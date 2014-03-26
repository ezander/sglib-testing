%function rv_map_approx
clc
%%
rand_seed(935590374)

m0 = 2;
m = 2;
V_y=gpcbasis_create('H', 'm', m0, 'p', 2);
y_j_beta = gpc_rand_coeffs(V_y, m);

n = 2;
V_x=gpcbasis_create('H', 'm', m0, 'p', 3);
x_i_alpha = gpc_rand_coeffs(V_x, n);

p_phi = 4;
p_int = 9;

[phi_j_gamma, V_phi]=mmse_estimate_gpc(x_i_alpha, V_x, y_j_beta, V_y, p_phi, p_int, 'cond_warning', 1e10);

%% show difference between approximations
xi = gpcgerm_sample(V_y);
y = gpc_evaluate(y_j_beta, V_y, xi);
x = gpc_evaluate(x_i_alpha, V_x, xi);

x_approx = gpc_evaluate(phi_j_gamma, V_phi, y);
underline('difference for one sample');
fprintf('y_%d: %6.4f  approx: %6.4f diff: %6.4f \n', [(1:n)', x, x_approx, abs(x-x_approx)]');

%%
N=100;
xi = gpcgerm_sample(V_y,N);
y = gpc_evaluate(y_j_beta, V_y, xi);
x = gpc_evaluate(x_i_alpha, V_x, xi);
x_approx = gpc_evaluate(phi_j_gamma, V_phi, y);

subplot(2,1,1)
%plot(x(1,:), x(2,:), 'b.', x_approx(1,:), x_approx(2,:), 'r.');
plot(x(1,:), x(2,:), 'ko', x_approx(1,:), x_approx(2,:), 'r.');

%% approximate L2 error of approximation by MC (LHS) sampling
subplot(2,1,2)
N=1000;
[y_relerr, xi] = compute_mc_error(x_i_alpha, V_x, y_j_beta, V_y, phi_j_gamma, V_phi, N);
plot(xi(1,:), xi(2,:), '.')
underline('L2 error');
fprintf('relerr y_%d: %6.4f%%\n', [(1:n)', y_relerr*100]');

