function Xn_i_alpha=mmse_update_gpc_basic(X_i_alpha, Y_func, V, ym, p_phi, p_int_mmse, p_int_proj)
% MMSE_UPDATE_GPC_BASIC Short description of mmse_update_gpc_basic.
%   MMSE_UPDATE_GPC_BASIC Long description of mmse_update_gpc_basic.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example mmse_update_gpc_basic">run</a>)
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

% Convert the GPC expansion of X into a function object
X_func = gpc_function(X_i_alpha, V);

% Now compute the MMSE estimator for X given Y and make a function
% out of this estimator
[phi_j_delta,V_phi]=mmse_estimate(X_func, Y_func, V, p_phi, p_int_mmse);
phi_func = gpc_function(phi_j_delta, V_phi);

% Create the prediction stochastic model for X as function
% composition between Y and phi and compute its GPC expansion
XM_func = funcompose(Y_func, phi_func);
XM_i_alpha = gpc_projection(XM_func, V, p_int_proj);

% Compute the best estimator value xm=phi(ym)
xm = funcall(phi_func, ym);

% Subtract the old coefficients to get the update and add xm
Xn_i_alpha =  X_i_alpha - XM_i_alpha;
Xn_i_alpha(:,1) = xm;

% The updated model xn and the update should be orthogonal
%assert(norm(gpc_covariance(xn_i_alpha, V_x, xm_i_alpha))<1e-10)
