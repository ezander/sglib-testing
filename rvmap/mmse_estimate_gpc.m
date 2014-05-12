function [phi_i_delta, V_phi]=mmse_estimate_gpc(X_i_alpha, V_X, Y_j_beta, V_Y, p_phi, p_int, varargin)
% MMSE_ESTIMATE_GPC Compute the MMSE estimator for GPC variables.
%   [PHI_I_DELTA, V_PHI]=MMSE_ESTIMATE_GPC(X_I_ALPHA, V_X, Y_J_BETA, V_Y,
%   P_PHI, P_INT, OPTIONS) computes the minimum mean square error estimator
%   PHI that minimises the error between X and PHI(Y). Here, X is given as
%   a GPC by X_I_ALPHA and V_X, Y by Y_J_BETA and V_Y. Both must be defined
%   on the same GPC germ. PHI is represented by multivariate polynomials of
%   maximum degree P_PHI. The coefficients are returned in PHI_I_DELTA, and
%   the system of polynomials if described by V_PHI (same for other GPC
%   functions). The one dimensional basis polynomials are the monomials by
%   default, but that can be changed using the POLYSYS option. P_INT is the
%   order of integration used.
%
% Options
%    'polysys': {'M'}, 'p', 'P', 'U', 'T', 'H', ...
%       The polynomial system used for representing PHI. Any sort of GPC
%       basis polynomial can be used.
%    'cond_warning': double, {inf}
%       A treshold value for the condition number of the linear system that
%       needs to be solved. If the condition number estimate is higher a
%       warning message is shown.
%
% References
%    [1] http://en.wikipedia.org/wiki/Minimum_mean_square_error
%    [2] D. P. Bertsekas and J. N. Tsitsiklis, Introduction to probability,
%        2 ed., Athena Scientific, 2008. 
%   MMSE_ESTIMATE_GPC Long description of mmse_estimate_gpc.
%
% Example (<a href="matlab:run_example mmse_estimate_gpc">run</a>)
%
% See also GPC, MMSE_ESTIMATE

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


% Check that bases of x and y are compatible
assert(gpcbasis_size(V_X, 2)==gpcbasis_size(V_Y, 2));
% check_gpc_compatibility(V_x, V_y, 'same_germ');

% Generate functions from the GPC bases and coeffs
X_func = gpc_function(X_i_alpha, V_X);
Y_func = gpc_function(Y_j_beta, V_Y);

% Call the MMSE estimator function
[phi_i_delta, V_phi]=mmse_estimate(X_func, Y_func, V_X, p_phi, p_int, varargin{:});
