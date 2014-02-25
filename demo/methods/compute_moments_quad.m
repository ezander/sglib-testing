function [u_mean, u_var, model] = compute_moments_quad(model, a_alpha, V_a, p_int, varargin)
% COMPUTE_MOMENTS_QUAD Compute mean and variance by high-dimensional quadrature.
%   [U_MEAN, U_VAR] = COMPUTE_MOMENTS_QUAD(INIT_FUNC, SOLVE_FUNC, P) computes
%   the mean and variance of a system described by INIT_FUNC and SOLVE_FUNC
%   by a quadrature rule of order P. The distribution is specified
%   by POLYSYS (see <a href="matlab:help gpc">GPC</a>).
%
% Options
%   grid: {'smolyak'}, 'full_tensor'
%     Choose the grid to use for integration points.
%
% Example (<a href="matlab:run_example compute_moments_quad">run</a>)
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

options = varargin2options(varargin);
[grid, options]=get_option(options, 'grid', 'smolyak');
check_unsupported_options(options, mfilename);

[x, w] = gpc_integrate([], V_a, p_int, 'grid', grid);
Q = length(w);

a = gpc_evaluate(a_alpha, V_a, x);
u = zeros(model.model_info.num_vars, Q);
for j = 1:Q
    a_j = a(:, j);
    [u(:, j), model] = model_solve(model, a_j);
end

u_mean = u * w;
u_var = binfun(@minus, u, u_mean).^2 * w;
