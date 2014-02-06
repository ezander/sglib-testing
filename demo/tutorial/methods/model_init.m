function state = model_init(num_params, num_vars, solve_func, step_func)
% MODEL_INIT Short description of model_init.
%   MODEL_INIT Long description of model_init.
%
% Options
%
% Example (<a href="matlab:run_example model_init">run</a>)
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


info = struct();
info.num_params = num_params;
info.num_vars = num_vars;
info.solve_func = solve_func;
info.step_func = step_func;

stats = struct();
stats.num_solve_calls = 0;
stats.num_step_calls = 0;
stats.time_solve_calls = 0;
stats.time_step_calls = 0;

state = struct();
state.model_info = info;
state.model_stats = stats;
