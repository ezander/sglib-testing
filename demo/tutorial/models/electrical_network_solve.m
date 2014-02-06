function [u, solve_info, state] = electrical_network_solve(state, p, varargin)
% ELECTRICAL_NETWORK_SOLVE Short description of electrical_network_solve.
%   ELECTRICAL_NETWORK_SOLVE Long description of electrical_network_solve.
%
% Example (<a href="matlab:run_example electrical_network_solve">run</a>)
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
[u0, options] = get_option(options, 'u0', zeros(state.model_info.num_vars,1));
check_unsupported_options(options, mfilename);

[u, iter, res, state] = general_iterative_solver(state, p, 'u0', u0, 'verbose', false);
solve_info.iter = iter;
solve_info.res = res;
