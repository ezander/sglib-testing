function plot_lmk3
% fig 2.3 page 25
clf
sig = 1;
N = 6;
[pos, els]=create_mesh_1d(0, 1, 40);
b=1;
[r_i_k, sigma_k] = kl_evp_analytic(pos, els, N, sig, b);
K = @(x,y)(sig^2 * exp(-abs(x-y)/b));
[x,y]=meshgrid(pos);

multiplot_init(3);

multiplot;
plot_3d(x,y,K(x,y));
title('$C_{U,U}$');

multiplot;
plot_3d(x,y,r_i_k*diag(sigma_k.^2)*r_i_k');
title('$C_{\hat{U},\hat{U}}$');

multiplot;
plot_3d(x,y,K(x,y)-r_i_k*diag(sigma_k.^2)*r_i_k');
title('$C_{U,U}-C_{\hat{U},\hat{U}}$');
    

function plot_3d(x,y,z)
surf(x,y,z);
hold all; 
contour(x,y,z);

