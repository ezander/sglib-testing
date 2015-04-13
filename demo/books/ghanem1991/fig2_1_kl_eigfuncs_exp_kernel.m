function fig2_1_kl_eigfuncs_exp_kernel
% Fig 2.1, page 29, (R. G. Ghanem & P. D. Spanos, Stochastic Finite
% Elements, 2nd edition from 2003)

clf

% Solve analytically
sig = 1;
N = 4;
b=1;
x = linspace(-0.5, 0.5);
[r_i_k, sigma_k] = kl_solve_1d_exp(x, sig, b, N);

clf
set(gca,'ColorOrder',[0,0,0])
set(gca,'LineStyleOrder',{'-', ':', '--', '-.'})
hold on
plot(x,r_i_k);
hold off
xlim([-0.5,1]);
ylim([-1.5,2.5]);
xlabel('LENGTH, x');
ylabel('EIGENFUNCTION,');
legend('First Eigenf.', 'Second Eigenf.', 'Third Eigenf.', 'Fourth Eigenf.');