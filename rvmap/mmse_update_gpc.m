function [pn_i_beta, V_pn]=mmse_update_gpc(p_func, y_func, V_py, ym, eps_func, V_eps, p_phi, p_int_mmse, p_pn, p_int_proj)
% MMSE_UPDATE_GPC Short description of mmse_update_gpc.
%   MMSE_UPDATE_GPC Long description of mmse_update_gpc.
%
% Options
%
% References
%
% Notes
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

% Combine the bases of the a priori model (V_py) and of the error model
% (V_eps)
[V_m, P_py, P_eps, xi_xy_ind, xi_eps_ind] = gpcbasis_combine(V_py, V_eps, 'outer_sum');

% Create evaluation function for a) measurement + error model (ym_func) and
% b) for the parameter to be estimated (pm_func), but now as functions in
% the combined space V_m
ym_func = @(xi)(funcall(y_func, xi(xi_xy_ind,:)) + funcall(eps_func, xi(xi_eps_ind, :)));
pm_func = @(xi)(funcall(p_func, xi(xi_xy_ind,:)));

% Now compute the MMSE estimator for P given Y+error and make a function
% out of this estimator
[phi_j_delta,V_phi]=mmse_estimate(pm_func, ym_func, V_m, p_phi, p_int_mmse);
phi_func = gpc_function(phi_j_delta, V_phi);

% Create the new stochastic model for the parameter as function composition
% between Y and phi and compute its GPC expansion
pn_func = funcompose(ym_func, phi_func);
V_pn = gpcbasis_create(V_m, 'p', p_pn);
pn_i_beta = gpc_projection(pn_func, V_pn, p_int_proj);

% Replace the mean value in the GPC of P by the best estimator value
% pm=phi(ym)
pm = funcall(phi_func, ym);
pn_i_beta(:,1) = pm;
