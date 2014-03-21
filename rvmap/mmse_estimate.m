function [phi_i_delta, V_phi]=mmse_estimate(x_i_alpha, V_x, y_j_beta, V_y, p_phi, p_int, varargin)
% APPROX_RVMAP Short description of approx_rvmap.
%   APPROX_RVMAP Long description of approx_rvmap.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example approx_rvmap">run</a>)
%
% See also

%   Elmar Zander
%   Copyright 2013, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

options=varargin2options(varargin);
[cond_warning,options]=get_option(options, 'cond_warning', false);
check_unsupported_options(options, mfilename);

% Check that bases of x and y are compatible
assert(gpcbasis_size(V_y, 2)==gpcbasis_size(V_x, 2));
% check_gpc_compatibility(V_x, V_y, 'same_germ');

% Determine dimension of y can create function basis
m = size(y_j_beta, 1);
V_phi=gpcbasis_create('U', 'm', m, 'p', p_phi);

% Generate integration points
[xi_k, w_k] = gpc_integrate([], V_y, p_int);

% Evaluate x, y and phi(y) at the integration points
x_i_k = gpc_evaluate(x_i_alpha, V_x, xi_k);
y_j_k = gpc_evaluate(y_j_beta, V_y, xi_k);
phi_gamma_k = gpcbasis_evaluate(V_phi, y_j_k);

% Compute matrix A and right hand side b and solve
phiw_gamma_k = binfun(@times, phi_gamma_k, w_k');
A = phiw_gamma_k * phi_gamma_k';
b = x_i_k * phiw_gamma_k';
phi_i_delta = (A\b)';
size(A)
if cond_warning
    kappa = condest(A);
    if kappa>=cond_warning
        warning('sglib:mmse_estimate:cond', ...
            'Condition number of matrix too large (%g), function approximation may be inaccurate', ...
            kappa);
    end
end
