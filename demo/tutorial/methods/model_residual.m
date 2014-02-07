function [res, minfo] = model_residual(minfo, curr_sol, params, varargin)
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

res_func = minfo.model_info.res_func;
%start = tic;
[res, minfo]=funcall(res_func, minfo, curr_sol, params, varargin{:});
%t = toc(start);
% minfo.model_stats.num_res_calls = minfo.model_stats.num_res_calls + 1;
% minfo.model_stats.time_res_calls = minfo.model_stats.time_res_calls + t;
%minfo.res_info = res_info;
