function [rmse, rmse_std, model] = compute_mc_error(model, a_alpha, V_a, u_alpha, V_u, N_or_xi, varargin)
% COMPUTE_MC_ERROR Short description of compute_mc_error.
%   COMPUTE_MC_ERROR Long description of compute_mc_error.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example compute_mc_error">run</a>)
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

options=varargin2options(varargin);
[mode,options]=get_option(options, 'mode', '');
[use_surrogate,options]=get_option(options, 'use_surrogate', true);
[solve_opts,options]=get_option(options, 'solve_opts', {});
[M,options]=get_option(options, 'repeat', 1);
check_unsupported_options(options, mfilename);

if isscalar(N_or_xi)
    sample_options=struct;
    if ~isempty(mode)
        sample_options.mode = mode;
    end
    N = N_or_xi;
    xi = gpcgerm_sample(V_a, N*M, sample_options);
else
    xi = N_or_xi;
end


a = gpc_evaluate(a_alpha, V_a, xi);

rmse_arr=nan(M,1);
for m=1:M
    mse = nan(N,1);
    for j = 1:N
        ind = j + (m-1)*N;
        a_j = a(:,ind);
        u_j_ap = gpc_evaluate(u_alpha, V_u, xi(:, ind));
        if use_surrogate
            [u_j_ex, model] = model_solve(model, a_j, 'u0', u_j_ap, solve_opts{:});
        else
            [u_j_ex, model] = model_solve(model, a_j, solve_opts{:});
        end
        mse(j) = sum((u_j_ex-u_j_ap).^2);
    end
    rmse_arr(m) = sqrt(mean(mse));
end
rmse = mean(rmse_arr);
rmse_std = std(rmse_arr);

% if nargout>1
%     mse = sort(mse);
%     n1 = round(N * 0.025);
%     n2 = round(N * 0.975);
%     rmse_int = reshape(sqrt(mse([n1, n2])), 1, []);
% end
