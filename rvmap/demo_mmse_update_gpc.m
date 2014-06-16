function demo_mmse_update_gpc
% DEMO_MMSE_UPDATE_GPC Show that the linear MMSE estimator is reproduced.
%
% References
%   [1] http://en.wikipedia.org/wiki/Minimum_mean_square_error#Linear_MMSE_estimator
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

linear_mmse_estimator_for_observation_process



function linear_mmse_estimator_for_observation_process

rand_seed(1234553412345)

% Generate some random GPC expansions for X and Z (mean Z must be zero!)
V_x = gpcbasis_create('H', 'm', 4, 'p', 1);
x_i_alpha = gpc_rand_coeffs(V_x, 4);

V_eps = gpcbasis_create('H', 'm', 5, 'p', 1);
eps_j_beta = gpc_rand_coeffs(V_eps, 3, 'zero_mean', true);
eps_func = gpc_function(eps_j_beta, V_eps);

% Model for Y is y=Ax+z (we could do that by function evaluation, but we
% can also directly combine the GPCs)
A = rand(3, 4);
y_func = gpc_function(A * x_i_alpha, V_x);

% Assume some concrete measurements for y in ym
ym = rand(3,1);

% Do the update via the MMSE method
p_phi=1;
p_int_mmse=2;
p_xn=1;
p_int_proj=2;
[xn_i_beta, V_xn]=mmse_update_gpc(x_i_alpha, y_func, V_x, ym, eps_func, V_eps, p_phi, p_int_mmse, p_xn, p_int_proj);


% Do the classical linear update (minimising y - (Wx+b))
Cx = gpc_covariance(x_i_alpha, V_x);
Cz = gpc_covariance(eps_j_beta, V_eps);

mean_x = gpc_moments(x_i_alpha, V_x);
mean_y = gpc_integrate(y_func, V_x, 2);

W=(Cx*A') / (A*Cx*A'+Cz);
b = mean_x - W * mean_y;

% This would be the classical new estimate (xn=W*ym+b) and covariance matrix Cxn
Cxn = W * (A*Cx);
xn = W * (ym - mean_y) + mean_x;
Ce = Cx - Cxn;


% Compare to the from 
assert_equals(Ce, gpc_covariance(xn_i_beta, V_xn), 'Ce match')
assert_equals(xn, gpc_moments(xn_i_beta, V_xn), 'xn match')
assert_equals(W * ym + b, gpc_moments(xn_i_beta, V_xn), 'xn match (2)')

