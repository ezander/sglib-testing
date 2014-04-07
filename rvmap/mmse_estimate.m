function [phi_i_delta, V_phi]=mmse_estimate(x_func, y_func, V_x, p_phi, p_int, varargin)
% MMSE_ESTIMATE Compute the MMSE estimator.
%   [PHI_I_DELTA, V_PHI]=MMSE_ESTIMATE(X_FUNC, Y_FUNC, V_X, P_PHI, P_INT,
%   OPTIONS) computes the minimum mean square error estimator PHI that
%   minimises the error between X and PHI(Y). Here, X is given by the
%   function X_FUNC, Y by the function Y_FUNC, and both are defined on the
%   same GPC germ V_X. PHI is represented by multivariate polynomials of
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
%
% Example (<a href="matlab:run_example mmse_estimate">run</a>)
%    
% See also GPC, MMSE_ESTIMATE_GPC

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

options=varargin2options(varargin);
[cond_warning,options]=get_option(options, 'cond_warning', inf);
[polysys,options]=get_option(options, 'polysys', 'M');
check_unsupported_options(options, mfilename);

% Generate integration points
[xi_k, w_k] = gpc_integrate([], V_x, p_int);

% Evaluate x and y at the integration points
x_i_k = funcall(x_func, xi_k);
y_j_k = funcall(y_func, xi_k);

% Determine dimension of codomain of y and create function basis
m = size(y_j_k, 1);
V_phi=gpcbasis_create(polysys, 'm', m, 'p', p_phi);
phi_gamma_k = gpcbasis_evaluate(V_phi, y_j_k);

% Compute matrix A and right hand side b and solve
phiw_gamma_k = binfun(@times, phi_gamma_k, w_k');
A = phiw_gamma_k * phi_gamma_k';
b = x_i_k * phiw_gamma_k';
phi_i_delta = (A\b')';

% Issue warning if the condition number is too high
if isfinite(cond_warning)
    kappa = condest(A);
    if kappa>=cond_warning
        warning('sglib:mmse_estimate_gpc:cond', ...
            'Condition number of matrix too large (%g), function approximation may be inaccurate', ...
            kappa);
    end
end
