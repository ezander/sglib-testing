function fig2_3_approx_correlation_func
% fig 2.3 page 25
clf
sig = 1;
N = 6;
[pos, ~]=create_mesh_1d(0, 1, 40);
b=1;
[r_i_k, sigma_k] = kl_solve_1d_exp(pos, sig, b, N);
K = @(x,y)(sig^2 * exp(-abs(x-y)/b));
[x,y]=meshgrid(pos);

R_uu = K(x,y);
Rk_uu = r_i_k*diag(sigma_k.^2)*r_i_k';

multiplot_init(3);
set(gcf, 'defaulttextinterpreter', 'latex');

multiplot;
plot_3d(x,y,R_uu);
title('$C_{U,U}$');

multiplot;
plot_3d(x,y,Rk_uu);
title('$C_{\hat{U},\hat{U}}$');

multiplot;
plot_3d(x,y,R_uu-Rk_uu);
title('$C_{U,U}-C_{\hat{U},\hat{U}}$');
    

function plot_3d(x,y,z)
surf(x,y,z);
hold on;
[~, hh] = contour3(x, y, z, 20);
        
% size zpos to match the data
zlims = get(gca, 'ZLim');
zpos = zlims(1);
for i = 1 : length(hh)
    zz = get(hh(i), 'ZData');
    set(hh(i), 'ZData', zpos * ones(size(zz)));
end

