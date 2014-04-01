function demo_mmse_linear
% DEMO_MMSE_LINEAR Show that the linear MMSE estimator is reproduced.
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

linear_mmse_estimator
linear_mmse_estimator_for_observation_process

function linear_mmse_estimator

% Generate some random GPC expansions for X and Y
V_x = gpcbasis_create('u', 'm', 4, 'p', 1);
x_i_alpha = gpc_rand_coeffs(V_x, 4);
V_y = V_x;
y_i_alpha = gpc_rand_coeffs(V_x, 3);


% Do the update via the MMSE method
[phi_i_gamma, V_phi] = mmse_estimate_gpc(x_i_alpha, V_x, y_i_alpha, V_y, ...
    1, 2, 'polysys', 'M');

% Do the classical linear update (minimising y - (Wx+b))
Cyx = gpc_covariance(y_i_alpha, V_x, x_i_alpha);
Cyy = gpc_covariance(y_i_alpha, V_y);
mean_x = gpc_moments(x_i_alpha, V_x);
mean_y = gpc_moments(y_i_alpha, V_y);

W = (Cyy\Cyx)';
b = mean_x - W * mean_y;

% Compare to the MMSE result (Note: do a comparison like this it is
% necessary that the ordering of the functions in V_phi is like 1, x_1,
% x_2, ..., x_m; otherwise the functions someselves would have to compared.
% To make sure this is the case polysys='M' has been chosen in the call to
% mmse_estimate_gpc).
assert_equals(phi_i_gamma, [b, W], 'phi(x)=Wx+b')



function linear_mmse_estimator_for_observation_process

% Generate some random GPC expansions for X and Z (mean Z must be zero!)
V_x = gpcbasis_create('H', 'm', 4, 'p', 1);
x_i_alpha = gpc_rand_coeffs(V_x, 4);

V_z = gpcbasis_create('H', 'm', 5, 'p', 1);
z_j_beta = gpc_rand_coeffs(V_z, 3);
z_j_beta(:,1) = 0;

% Combine stuff into one probability space
[xz_ij_alpha, V_xz] = gpc_combine_inputs(x_i_alpha, V_x, z_j_beta, V_z);
x_i_alpha = xz_ij_alpha(1:4, :);
z_j_alpha = xz_ij_alpha(5:end, :);

% Covariance between X and Z should be zero
Cxz = gpc_covariance(x_i_alpha, V_xz, z_j_alpha);
assert_equals(Cxz, zeros(4,3), 'Cxy=0')

% Model for Y is y=Ax+z (we could do that by function evaluation, but we
% can also directly combine the GPCs)
A = rand(3, 4);
y_j_alpha = A * x_i_alpha + z_j_alpha;

% Do the update via the MMSE method
[phi_i_gamma, V_phi] = mmse_estimate_gpc(x_i_alpha, V_xz, y_j_alpha, V_xz, ...
    1, 2, 'polysys', 'M');

% Do the classical linear update (minimising y - (Wx+b))
Cxx = gpc_covariance(x_i_alpha, V_xz);
Cyx = gpc_covariance(y_j_alpha, V_xz, x_i_alpha);
Cyy = gpc_covariance(y_j_alpha, V_xz);
Czz = gpc_covariance(z_j_alpha, V_xz);
assert_equals(Cyy, A*Cxx*A'+Czz, 'Cyy=A*Cxx*A''+Czz');
assert_equals(Cyx, A*Cxx, 'Cyx=A*Cxx');

mean_x = gpc_moments(x_i_alpha, V_xz);
mean_y = gpc_moments(y_j_alpha, V_xz);

W=(Cxx*A') * inv(A*Cxx*A'+Czz);
b = mean_x - W * mean_y;

% Compare to the MMSE result 
assert_equals(phi_i_gamma, [b, W], 'phi(x)=Wx+b')
