function [u_i_alpha, state, x, w] = compute_response_surface_projection(state, a_alpha, V_a, V_u, p_int, varargin)
% COMPUTE_RESPONSE_SURFACE_PROJECTION Compute a gpc response surface representation.
%   COMPUTE_RESPONSE_SURFACE_PROJECTION Long description of compute_response_surf_projection.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example compute_response_surface_projection">run</a>)
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
[grid,options]=get_option(options, 'grid', 'smolyak');
check_unsupported_options(options, mfilename);

[x,w] = gpc_integrate([], V_u, p_int, 'grid', grid);
M = gpcbasis_size(V_u, 1);
Q = length(w);
u_i_alpha = zeros(state.model_info.num_vars, M);
a = gpc_evaluate(a_alpha, V_a, x);

for j = 1:Q
    a_j = a(:, j);
    [u_i_j, state] = model_solve(state, a_j);
    x_j = x(:, j);
    psi_j_alpha_dual = gpcbasis_evaluate(V_u, x_j, 'dual', true);
    u_i_alpha = u_i_alpha + w(j) * u_i_j * psi_j_alpha_dual;
end
