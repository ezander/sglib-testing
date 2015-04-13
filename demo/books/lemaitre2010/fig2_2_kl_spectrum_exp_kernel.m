function fig2_2_kl_spectrum_exp_kernel
% fig 2.2 page 23

clf
N = 20;
sig = 1;
[pos, ~]=create_mesh_1d(0, 1, 2);
hold off;
multiplot_init(1,2);

multiplot;
for b=linspace(0.1, 1, 10)
    [~, sigma_k] = kl_solve_1d_exp(pos, sig, b, N);
    loglog(1:N, sigma_k.^2, '-'); hold all;
    xlim([1,N])
end

multiplot;
for b=linspace(1, 10, 10)
    [~, sigma_k] = kl_solve_1d_exp(pos, sig, b, N);
    loglog(1:N, sigma_k.^2, '-'); hold all;
    xlim([1,N])
end
hold off



