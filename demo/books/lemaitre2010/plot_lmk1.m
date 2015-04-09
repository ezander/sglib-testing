function plot_lmk1
% fig 2.1 page 22 (lemaitre & knio spectral methods)

clf
multiplot_init(2,2,'ordering','row');

sig = 1;
N = 10;
b=1;
[pos, els]=create_mesh_1d(0, 1, 100);
[r_i_k, sigma_k] = kl_evp_analytic(pos, els, N, sig, b);


multiplot; 
plot(pos, binfun(@times, r_i_k, sigma_k));

multiplot; 
semilogy(1:N, sigma_k.^2,'x');
grid on;

[r_i_k, sigma_k] = kl_evp_numerically(pos, els, N, sig, b);
r_i_k = binfun(@times, r_i_k, sign(r_i_k(52,:)));

multiplot; 
plot(pos, binfun(@times, r_i_k, sigma_k));

multiplot; 
semilogy(1:N, sigma_k.^2,'x');
grid on;




