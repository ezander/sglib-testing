function fig2_1_kl_exp_kernel
% fig 2.1 page 22 (lemaitre & knio spectral methods)

clf
multiplot_init(2, 2, 'ordering', 'row');

% Solve analytically
sig = 1;
N = 10;
b=1;
[pos, els]=create_mesh_1d(0, 1, 100);
[r_i_k, sigma_k] = kl_solve_1d_exp(pos, 1, b, N);

multiplot; 
plot(pos, binfun(@times, r_i_k, sigma_k));

multiplot; 
semilogy(1:N, sigma_k.^2,'x');
grid on;

% Solve numerically (second line corrects for differences in sign of the
% functions)
[r_i_k, sigma_k] = kl_evp_numerically(pos, els, N, sig, b);
r_i_k = binfun(@times, r_i_k, sign(r_i_k(52,:)));

multiplot; 
plot(pos, binfun(@times, r_i_k, sigma_k));

multiplot; 
semilogy(1:N, sigma_k.^2,'x');
grid on;




function [r_i_k, sigma_k] = kl_evp_numerically(pos, els, N, sig, b)
G_N=mass_matrix(pos,els);
C=covariance_matrix(pos,...
    funcreate(@exponential_covariance, @funarg, @funarg, b, sig));
[r_i_k,sigma_k]=kl_solve_evp( C, G_N, N);

