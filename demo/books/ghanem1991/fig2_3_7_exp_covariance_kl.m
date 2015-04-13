function fig2_3_7_exp_covariance_kl
% Fig 2.3-2.7, pages 30-32, (R. G. Ghanem & P. D. Spanos, Stochastic Finite
% Elements, 2nd edition from 2003)
% fig 2.3 page 25
clf
sig = 1;
a = 0.5;
b = 1;
x=linspace(-a,a,30);
[X,Y]=meshgrid(x);

% Exact covariance
K = @(x,y)(sig^2 * exp(-abs(x-y)/b));
R_uu = K(X,Y);

multiplot_init(3, 2, 'ordering', 'row');
set(gcf, 'defaulttextinterpreter', 'latex');

multiplot;
plot_3d(X,Y,R_uu);
multiplot;

% 4 term approximation
N = 4;
[r_i_k, sigma_k] = kl_solve_1d_exp(x, sig, b, N);
Rk_uu = r_i_k*diag(sigma_k.^2)*r_i_k';

multiplot;
plot_3d(X,Y,Rk_uu);

multiplot;
plot_3d(X,Y,R_uu-Rk_uu);

% 10 term approximation
N = 10;
[r_i_k, sigma_k] = kl_solve_1d_exp(x, sig, b, N);
Rk_uu = r_i_k*diag(sigma_k.^2)*r_i_k';

multiplot;
plot_3d(X,Y,Rk_uu);

multiplot;
plot_3d(X,Y,R_uu-Rk_uu);


function plot_3d(x,y,z)
mesh(x,y,z,'EdgeColor','k');
