%% Represent a Gaussian random field
clear
clf

[cov_func, x, els, l_c] = setup_model();
N = size(x, 2);

% Compute the KL expansion
L = min(30, N-1);
C = covariance_matrix(x, cov_func);
G = mass_matrix(x, els);
[r_i_k, sigma_k] = kl_solve_evp(C, G, L);

% Setup the multiindex for m = L Gaussian random variables
m = L;
I = multiindex(m, 1);

% The corresponding coefficients have 
k = 1:m;
theta_k_alpha = zeros(L,m+1);
theta_k_alpha(k,k+1) = eye(m);

% Convert to compact form and include mean
% [u2_i_k, theta2_k_alpha] = kl_pce_to_compact_form(sin(pi*x'), r_i_k, sigma_k, theta_k_alpha);

% The conversion made explicit
% a) Multiply in the sigmas 
u_i_k = binfun(@times, r_i_k, sigma_k);

% b) Put the mean inside
u_i_k = [sin(pi*x'), u_i_k];
theta_k_alpha = [1, zeros(1,m); theta_k_alpha];


% Generate some normally distributed random numbers for plotting some
% realisations of the random field
xi = randn(m, 40);
clf
plot(x, kl_pce_field_realization(u_i_k, theta_k_alpha, I, xi))

% Compute the mean and standard deviation of the KL and add it to the plot
[mu_r, sigma_r] = kl_pce_moments(u_i_k, theta_k_alpha, I);
hold all;
plot(x, mu_r, 'b', 'LineWidth', 2);
plot(x, [mu_r + sigma_r, mu_r - sigma_r], 'b-.', 'LineWidth', 2);
hold off;

%plot_kl_pce_realizations_1d( x, u_i_k, theta_k_alpha, I, 'colormap', 'cool', 'realizations', 100);
