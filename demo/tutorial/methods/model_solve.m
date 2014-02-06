function [sol, state] = model_solve(state, params, varargin)
% MODEL_SOLVE Short description of model_solve.
%   MODEL_SOLVE Long description of model_solve.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example model_solve">run</a>)
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

solve_func = state.model_info.solve_func;
start = tic;
[sol, solve_info, state]=funcall(solve_func, state, params, varargin);
t = toc(start);
state.model_stats.num_solve_calls = state.model_stats.num_solve_calls + 1;
state.model_stats.time_solve_calls = state.model_stats.time_solve_calls + t;
state.solve_info = solve_info;
