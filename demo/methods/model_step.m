function [sol, model] = model_step(model, old_sol, params, varargin)
% MODEL_STEP Short description of model_step.
%   MODEL_STEP Long description of model_step.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example model_step">run</a>)
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


step_func = model.model_info.step_func;
start = tic;
[sol, model]=funcall(step_func, model, old_sol, params, varargin{:});
t = toc(start);
model.model_stats.num_step_calls = model.model_stats.num_step_calls + 1;
model.model_stats.time_step_calls = model.model_stats.time_step_calls + t;
