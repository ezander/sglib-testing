function [u_i_alpha, model, x, w]=compute_response_surface_nonintrusive_galerkin(model, a_alpha, V_a, V_u, p_int, varargin)
% COMPUTE_RESPONSE_SURFACE_NONINTRUSIVE_GALERKIN Short description of compute_response_surface_nonintrusive_galerkin.
%   COMPUTE_RESPONSE_SURFACE_NONINTRUSIVE_GALERKIN Long description of compute_response_surface_nonintrusive_galerkin.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example compute_response_surface_nonintrusive_galerkin">run</a>)
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
[u0_i_alpha, options]=get_option(options, 'initial_sol', []);
[maxiter, options]=get_option(options, 'maxiter', 50);
[steptol, options]=get_option(options, 'steptol', 1e-5);
[grid, options]=get_option(options, 'grid', 'smolyak');
[method, options]=get_option(options, 'method', 'simple');
[verbosity, options]=get_option(options, 'verbosity', 0);
check_unsupported_options(options, mfilename);

if isempty(u0_i_alpha)
    u0_i_alpha = zeros(model.model_info.num_vars, gpcbasis_size(V_u, 1));
    % TODO correct this
    % u_i_alpha = repmat(model.u0, 1, gpcbasis_size(V_u,1));
end

[x,w] = gpc_integrate([], V_u, p_int, 'grid', grid);
a=gpc_evaluate(a_alpha, V_a, x);

switch method
    case 'simple'
        [u_i_alpha, model]=nisg_simple_iterations(model, u0_i_alpha, V_u, a, x, w, maxiter, steptol, verbosity);
    case 'bfgs'
        [u_i_alpha]=nisg_bfgs_solve(u0_i_alpha, V_u, a, x, w, maxiter, steptol, verbosity);
end


    function u_i_alpha=nisg_bfgs_solve(u0_i_alpha, V_u, a, x, w, maxiter, steptol, verbosity)
        
        % solve by quasi newton stuff
        %func = funcreate(@nisg_step_bfgs, model, @funarg, V_u, a, x, w);
        func = funcreate(@nisg_step_bfgs, @funarg, V_u, a, x, w);
        
        qnewton_opts.abstol = steptol;
        qnewton_opts.verbosity = 1;
        qnewton_opts.maxiter= 100;
        qnewton_opts.line_search_func = @line_search_gradient;
        qnewton_opts.line_search_opts = {'alpha0', 1};
        %[func, stats] = CallStatistics.wrap_function(func);
        I=IdentityOperator.from_vector(u0_i_alpha);
        H0 = inv(LBFGSOperator(I, 5));
        [u_i_alpha, flag, info] = minfind_quasi_newton(func, u0_i_alpha, H0, qnewton_opts);
        % flag
        % info
    end


    function [u_i_alpha, model, x, w]=nisg_simple_iterations(model, u_i_alpha, V_u, a, x, w, maxiter, steptol, verbosity);
        converged=false;
        for k=1:maxiter
            % do the computation of unext
            [unext_i_alpha, model] = nisg_step(model, u_i_alpha, V_u, a, x, w);
            
            % compute the size of the update step in the Frobenius norm, maybe not
            % the best convergence criterion, but at least easy and fast to compute
            diff=norm(u_i_alpha-unext_i_alpha,'fro');
            u_i_alpha = unext_i_alpha;
            
            % print
            if verbosity>0
                strvarexpand('iter: $k$, diff: $diff$');
            end
            
            % check for convergence
            if diff<steptol
                converged = true;
                break;
            end
            
            %
            if ~isfinite(diff) || isnan(diff)
                break
            end
        end
        
        if ~converged
            warning('sglib:non_intr_galerkin', 'The Galerkin iteration did not converge after %d iterations.', maxiter);
        end
    end

    function [not_computed, unext_i_alpha] = nisg_step_bfgs(u_i_alpha, V_u, a, x, w)
        [unext_i_alpha, model] = nisg_step(model, u_i_alpha, V_u, a, x, w);
        unext_i_alpha = tensor_add(u_i_alpha, unext_i_alpha, -1);
        not_computed=nan;
        disp('call nisg_step_bfgs');
    end

    function [unext_i_alpha, model] = nisg_step(model, u_i_alpha, V_u, a, x, w)
        unext_i_alpha = u_i_alpha;
        % do the computation of unext
        
        % that's the z-sum here, i.e. z=x(:,j)
        Q = length(w);
        for j=1:Q
            x_j = x(:, j);
            a_j = a(:, j);
            % compute u_i_p = sum u_i_alpha Psi_alpha(p)
            u_i_p = gpc_evaluate(u_i_alpha, V_u, x_j);
            % evaluate S at p, u_i_p
            [S_p, model] = model_step(model, u_i_p, a_j);
            S_p = S_p - u_i_p;
            % update unext
            unext_i_alpha = unext_i_alpha + w(j) * S_p * gpcbasis_evaluate(V_u, x_j, 'dual', true);
        end
    end
end