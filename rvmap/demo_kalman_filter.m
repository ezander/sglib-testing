function demo_kalman_filter
% DEMO_KALMAN_FILTER Show how to reproduce the Kalman filter.
%
%   This demos tries to directly reproduce the results from [1], first in
%   the classical "Kalman manner", and then via a GPC model and the MMSE
%   functions. The naming of the variables closely follows [1].
%
% References
%   [1] http://en.wikipedia.org/wiki/Kalman_filter
%
% Example (<a href="matlab:run_example demo_kalman_filter">run</a>)
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

% F is the state transition model which is applied to the previous state x
% B is the control-input model which is applied to the control vector u
% w is the process noise which is assumed to be drawn from a zero mean
%   multivariate normal distribution with covariance Q. 
% H is the observation model which maps the true state space into the
%   observed space 
% v is the observation noise which is assumed to be zero mean Gaussian
%   white noise with covariance R

% State transition model:
%   xn = F * x + B * u + w
% Observation model:
%   z = H * x + v

%% Problem setup

rand_seed(1234553412345)
% Define the dimensions of the state x, control-input u, and observations z
nx = 5;
nu = 2;
nz = 3;

% Define the number of random variables x, w, and v are modelled upon
mx = 6;
mw = 4;
mv = 4;

% Now define some random models for them
V_x = gpcbasis_create('H', 'm', mx, 'p', 1);
x_i_alpha = gpc_rand_coeffs(V_x, nx);

V_w = gpcbasis_create('H', 'm', mw, 'p', 1);
w_i_beta = gpc_rand_coeffs(V_w, nx, 'zero_mean', true);

V_v = gpcbasis_create('H', 'm', mv, 'p', 1);
v_j_gamma = gpc_rand_coeffs(V_v, nz, 'zero_mean', true);

% Recover the covariance matrices from here
P = gpc_covariance(x_i_alpha, V_x);
Q = gpc_covariance(w_i_beta, V_w);
R = gpc_covariance(v_j_gamma, V_v);

% Also define some random matrices for the other stuff (i.e. the transition
% matrix F, the control matrix B, the observation matrix H, the control
% variable u, the actual measurement z) 
F = rand(nx, nx);
B = rand(nx, nu);
H = rand(nz, nx);
u = rand(nu, 1);
z = rand(nz, 1);
I = eye(nx);

%% Now the classical Kalman stuff
%  a) Prediction step (predicted new state xnp and covariance Pnp)

mean_x = gpc_moments(x_i_alpha, V_x);
xnp = F * mean_x + B * u;
Pnp = F * P * F' + Q;

%  b) Update step (innovation y, innovation covariance S, Kalman gain K,
%  updated state estimate xn, updated covariance Pn)
y = z - H * xnp;
S = H * Pnp * H' + R;
K = Pnp * H' / S;
xn = xnp + K * y;
Pn = (I - K*H)*Pnp;

%% And now with the MMSE functions on GPC variables

% Define a priori model for xn: 
%   xn = F * x + B * u + w
% (The tricky thing is that the separate probalility spaces for x and w
% have to be combined into one, and during the function evaluation split up
% again, so that the x and w part only get "their" variables)
[V_xnp, ~, ~, xi_x_ind, xi_w_ind] = gpcbasis_combine(V_x, V_w, 'outer_sum');
xnp_func = @(xi)(...
    F * gpc_evaluate(x_i_alpha, V_x, xi(xi_x_ind, :)) + ...  % F * x +
    repmat(B*u, 1, size(xi,2)) + ...                         % B * u +
    gpc_evaluate(w_i_beta, V_w, xi(xi_w_ind, :)));           % w

% Just another way to do it, by combining the GPC expansions and then
% making a function out of it
[xw_ii_gamma, V_xnp] = gpc_combine_inputs(x_i_alpha, V_x, w_i_beta, V_w);
xnp_i_gamma = [F, I] * xw_ii_gamma;
xnp_i_gamma(:,1) = xnp_i_gamma(:,1) + B*u;
xnp_func = gpc_function(xnp_i_gamma, V_xnp);

% In-between check that the predicted state xnp and covariance Pnp match
assert_equals(gpc_moments(xnp_i_gamma, V_xnp), xnp, 'predicted state estimate xnp')
assert_equals(gpc_covariance(xnp_i_gamma, V_xnp), Pnp, 'predicted covariance estimate Pnp')

% Define observation model (without the v, that comes extra)
%   z = H * xn
y_func = @(xi)(H * funcall(xnp_func, xi));

% In-between check of the covariance matrices of the observation
assert_equals(gpc_covariance(H*xnp_i_gamma, V_xnp)+R, S, 'innovation covariance S')

% Define error model (must be defined separately from observation model)
v_func = gpc_function(v_j_gamma, V_v);

% Now call the MMSE update procedure
p_phi=1;
p_int_mmse=2;
p_xn=1;
p_int_proj=2;
[xn_i_beta, V_xn]=mmse_update_gpc(xnp_func, y_func, V_xnp, z, v_func, V_v, p_phi, p_int_mmse, p_xn, p_int_proj);

% Final check of Kalman estimate and covariance update
assert_equals( gpc_moments(xn_i_beta, V_xn), xn, 'updated state estimate xn');
assert_equals( gpc_covariance(xn_i_beta, V_xn), Pn, 'updated covariance estimate Pn');

% ... and you see that the classical Kalman filter is exactly reproduced.
% The MMSE is completely independent of degree and in the limit case of all
% degrees equal to 1 reproduces the Kalman filter (or linear observation
% processes). 
