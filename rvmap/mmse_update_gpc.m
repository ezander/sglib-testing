function [Xn_i_beta, V_Xn]=mmse_update_gpc(X_i_alpha, Y_func, V_X, ym, eps_func, V_eps, p_phi, p_int_mmse, p_pn, p_int_proj)
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

% Combine the bases of the a priori model (V_X) and of the error model
% (V_eps)
[V_Xn, ind_V_X, ~, ind_xi_X, ind_xi_eps] = gpcbasis_combine(V_X, V_eps, 'outer_sum');

% Create evaluation function for measurement + error model (YM_func) 
YM_func = @(xi)(funcall(Y_func, xi(ind_xi_X,:)) + funcall(eps_func, xi(ind_xi_eps, :)));

% Extend the GPC coefficients X_i_alpha to the larger basis V_Xn
X_i_beta = zeros(size(X_i_alpha,1), gpcbasis_size(V_Xn,1));
X_i_beta(:,ind_V_X) = X_i_alpha;

Xn_i_beta=mmse_update_gpc_basic(X_i_beta, YM_func, V_Xn, ym, p_phi, p_int_mmse, p_int_proj);


