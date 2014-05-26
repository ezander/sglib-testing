function [Xn_i_beta, V_Xn]=mmse_update_gpc(X_func, Y_func, V_X, ym, eps_func, V_eps, p_phi, p_int_mmse, p_pn, p_int_proj)
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
[V_XM, ~, ~, xi_x_ind, xi_eps_ind] = gpcbasis_combine(V_X, V_eps, 'outer_sum');

% Create evaluation function for a) measurement + error model (ym_func) and
% b) for the parameter to be estimated (x2_func), but now as functions in
% the combined space V_m
YM_func = @(xi)(funcall(Y_func, xi(xi_x_ind,:)) + funcall(eps_func, xi(xi_eps_ind, :)));
X2_func = @(xi)(funcall(X_func, xi(xi_x_ind,:)));

% Now compute the MMSE estimator for X given Y+eps and make a function
% out of this estimator
[phi_j_delta,V_phi]=mmse_estimate(X2_func, YM_func, V_XM, p_phi, p_int_mmse);
phi_func = gpc_function(phi_j_delta, V_phi);

% Create the prediction stochastic model for X as function
% composition between Y and phi and compute its GPC expansion
XM_func = funcompose(YM_func, phi_func);
V_Xn = gpcbasis_create(V_XM, 'p', p_pn);
XM_i_beta = gpc_projection(XM_func, V_Xn, p_int_proj);

% Subtract the old coefficients to get the update
x_i_beta = gpc_projection(X2_func, V_Xn, p_int_proj);
Xn_i_beta =  x_i_beta - XM_i_beta;

% Replace the mean value in the GPC of XN by the best estimator value
% xm=phi(ym)
xm = funcall(phi_func, ym);
Xn_i_beta(:,1) = xm;

% The new model pn and the update should be orthogonal
CXn = norm(gpc_covariance(Xn_i_beta, V_Xn), 'fro');
CXnXM = norm(gpc_covariance(Xn_i_beta, V_Xn, XM_i_beta), 'fro');
if CXnXM>1e-10*CXn
    warning('sglib:mmse_update_gpc', 'gpc update not orthogonal (%g>1e-10*%g)', CXnXM, CXn);
end
