function demo_mmse_nonlinear_hgm
% DEMO_MMSE_NONLINEAR_HGM Short description of demo_mmse_nonlinear_hgm.
% Example (<a href="matlab:run_example demo_mmse_nonlinear_hgm">run</a>)
%
% See also

%   Elmar Zander
%   Copyright 2014, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.


alpha = 0.1;
delta = 0.04;
k = pi/2;
xi0 = pi/8;

x_func = @(xi)(cos(k*xi(1,:)));
y_func = @(xi)(k*xi(1,:) + alpha * xi(2,:));

V_x = gpcbasis_create('hh', 'p', 1);

p_phi=16; p_int=13;
[phi_i_delta, V_phi]=mmse_estimate(x_func, y_func, V_x, p_phi, p_int)
phi_func = gpc_function(phi_i_delta, V_phi);


xi = gpcgerm_sample(V_x, 3000);
x = funcall(x_func, xi);
y = funcall(y_func, xi);
x2 = funcall(phi_func,y);

ym = xi0 + delta;
xm = funcall(phi_func, ym)

mh=multiplot_init(2,2);

multiplot;
plot(x,y,'.', x2,y,'.',xm,ym,'rx',cos(k*xi0),(k*xi0), 'ro');
legend('a priori', 'pred', 'meas', 'truth');

multiplot;
plot_density(x);
hold on; plot(xm, 0, 'rx', 'markersize', 12); plot(cos(k*xi0), 0, 'ro', 'markersize', 12); hold off; 


xm_func = funcompose(y_func, phi_func);
xn1_func = @(xi)(xm + funcall(x_func, xi) - funcall(xm_func, xi));
xn1 = funcall(xn1_func, xi);

multiplot;
plot_density(xn1); 
hold on; plot(xm, 0, 'rx', 'markersize', 12); plot(cos(k*xi0), 0, 'ro', 'markersize', 12); hold off; 


xn2_func = @(xi)(funcall(x_func, xi) + ...
    funcall(phi_func, ym - funcall(y_func, xi)));
xn2 = funcall(xn2_func, xi);

multiplot;
plot_density(xn2);
hold on; plot(xm, 0, 'rx', 'markersize', 12); plot(cos(k*xi0), 0, 'ro', 'markersize', 12); hold off; 

same_scaling(mh, 'x');

save_figure(gcf, 'nonlinear_update', 'figdir', '.', 'fontsize', 12)

