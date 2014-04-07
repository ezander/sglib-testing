function [xn_i_beta, V_xn]=mmse_update_gpc(x_func, y_func, V_x, ym, eps_func, V_eps, p_phi, p_int_mmse, p_pn, p_int_proj)
% MMSE_UPDATE_GPC Update a GPC given some measurements and a measurement model.
%
% Example (<a href="matlab:run_example mmse_update_gpc">run</a>)
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

% Combine the bases of the a priori model (V_x) and of the error model
% (V_eps)
[V_xm, ~, ~, xi_x_ind, xi_eps_ind] = gpcbasis_combine(V_x, V_eps, 'outer_sum');

% Create evaluation function for a) measurement + error model (ym_func) and
% b) for the parameter to be estimated (x2_func), but now as functions in
% the combined space V_m
ym_func = @(xi)(funcall(y_func, xi(xi_x_ind,:)) + funcall(eps_func, xi(xi_eps_ind, :)));
x2_func = @(xi)(funcall(x_func, xi(xi_x_ind,:)));

% Now compute the MMSE estimator for X given Y+eps and make a function
% out of this estimator
[phi_j_delta,V_phi]=mmse_estimate(x2_func, ym_func, V_xm, p_phi, p_int_mmse);
phi_func = gpc_function(phi_j_delta, V_phi);

% Create the prediction stochastic model for X as function
% composition between Y and phi and compute its GPC expansion
xm_func = funcompose(ym_func, phi_func);
V_xn = gpcbasis_create(V_xm, 'p', p_pn);
xm_i_beta = gpc_projection(xm_func, V_xn, p_int_proj);

% Subtract the old coefficients to get the update
x_i_beta = gpc_projection(x2_func, V_xn, p_int_proj);
xn_i_beta =  x_i_beta - xm_i_beta;

% Replace the mean value in the GPC of XN by the best estimator value
% xm=phi(ym)
xm = funcall(phi_func, ym);
xn_i_beta(:,1) = xm;

% The new model pn and the update should be orthogonal
assert(norm(gpc_covariance(xn_i_beta, V_xn, xm_i_beta))<1e-10)
