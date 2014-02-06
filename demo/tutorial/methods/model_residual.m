function [res, state] = model_residual(state, curr_sol, params, varargin)
% MODEL_RESIDUAL Short description of model_residual.
%   MODEL_RESIDUAL Long description of model_residual.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example model_residual">run</a>)
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

res_func = state.model_info.res_func;
%start = tic;
[res, state]=funcall(res_func, state, curr_sol, params, varargin{:});
%t = toc(start);
% state.model_stats.num_res_calls = state.model_stats.num_res_calls + 1;
% state.model_stats.time_res_calls = state.model_stats.time_res_calls + t;
%state.res_info = res_info;
