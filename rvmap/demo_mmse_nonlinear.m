function demo_mmse_nonlinear
% DEMO_MMSE_NONLINEAR Short description of demo_mmse_nonlinear_hgm.
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


% Problem definition
alpha = 0.1;
delta = 0.04;
k = pi/2;
xi0 = pi/8;

V_x = gpcbasis_create('hh', 'p', 10);
x_func = @(xi)(cos(k*xi(1,:)));
[x_i_alpha] = gpc_projection(x_func, V_x, 10);
x_func=gpc_function(x_i_alpha,V_x);

y_func = @(xi)(k*xi(1,:) + alpha * xi(2,:));

% Define the actual "measurment"
ym = xi0 + delta;

% Get some samples for later use
xi = gpcgerm_sample(V_x, 3000);

% Definition of some algorithm parameters
p_phi=16; p_int=13;

% Compute the MMSE estimator
[phi_i_delta, V_phi]=mmse_estimate(x_func, y_func, V_x, p_phi, p_int)
phi_func = gpc_function(phi_i_delta, V_phi);

% Prepare the plotting
mh=multiplot_init(2,2);

% Sample from x and y and show the mapping phi
x = funcall(x_func, xi(:,1:3000));
y = funcall(y_func, xi(:,1:3000));
x2 = funcall(phi_func,y);
xm = funcall(phi_func, ym);

multiplot; plot(x,y,'.', x2,y,'.',xm,ym,'rx',cos(k*xi0),(k*xi0), 'ro');
legend('a priori', 'pred', 'meas', 'truth');

% Show the a priori density of x
multiplot;
x = funcall(x_func, xi);
plot_density(x, 'rug', true);
hold on; plot(xm, 0, 'rx', 'markersize', 12); plot(cos(k*xi0), 0, 'ro', 'markersize', 12); hold off; 


xn_i_alpha = mmse_update_gpc_basic(x_i_alpha, y_func, V_x, ym, p_phi, p_int, p_int);
xn = gpc_evaluate(xn_i_alpha, V_x, xi);

multiplot;
plot_density(xn, 'rug', true); 
hold on; plot(xm, 0, 'rx', 'markersize', 12); plot(cos(k*xi0), 0, 'ro', 'markersize', 12); hold off; 


same_scaling(mh, 'x');

save_figure(gcf, 'nonlinear_update', 'figdir', '.', 'fontsize', 12)

