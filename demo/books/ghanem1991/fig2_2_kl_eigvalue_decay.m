function fig2_2_kl_eigvalue_decay
% Fig 2.2, page 30, (R. G. Ghanem & P. D. Spanos, Stochastic Finite
% Elements, 2nd edition from 2003)


clf;
hold all;
set(gca,'ColorOrder',[0,0,0])
set(gca,'LineStyleOrder',{'-', ':', '--', '-.'})
b_vals = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0];

for b = b_vals
    sig = 1;
    N = 10;
    x = linspace(-0.5, 0.5);
    [~, sigma_k] = kl_solve_1d_exp(x, sig, b, N);
    lambda_k = sigma_k.^2;
    plot(lambda_k);
end

legend_parametric('b=%g', b_vals);
