function plot_lmk2
% fig 2.2 page 23

clf
sig = 1;
N = 20;
[pos, els]=create_mesh_1d(0, 1, 2);
hold off;
multiplot_init(1,2);

multiplot;
for b=linspace(0.1, 1, 10)
    [~, sigma_k] = kl_evp_analytic(pos, els, N, sig, b);
    loglog(1:N, sigma_k.^2, '-'); hold all;
    xlim([1,N])
end

multiplot;
for b=linspace(1, 10, 10)
    [~, sigma_k] = kl_evp_analytic(pos, els, N, sig, b);
    loglog(1:N, sigma_k.^2, '-'); hold all;
    xlim([1,N])
end
hold off



