%%
% initialise the minfo of the electrical network
clc
minfo = electrical_network_init('R', 0.02, 'f0', 0.1);
% display the minfo
minfo

% get the number of parameters from the minfo and set to some random value
p = rand(minfo.model_info.num_params, 1);
% get a starting vector from the minfo 
u = funcall(minfo.model_info.sol_init_func, p);
% compute the residual for the starting vector and the given minfo
res = electrical_network_residual(minfo, u, p);
norm(res)

%%
% solve the system with a nonlinear solver based on Picard iterations
minfo = electrical_network_init('R', 0.02, 'f0', 0.1, 'newton', false);
u = funcall(minfo.model_info.sol_init_func, p);
max_iter = 100;
abstol = 1e-5;

for iter = 1:max_iter
    res = model_residual(minfo, u, p);
    fprintf( 'iter %d, norm=%g\n', iter, norm(res));
    if norm(res)<abstol
        break;
    elseif iter==max_iter
        error('solve:no_conv', ...
            'Could not reach convergence (abstol=%g) in %d iterations', ...
            abstol, max_iter);
    end
    
    u = model_step(minfo, u, p);
end

%%
% the same as before, but solver has been put into the function
% 'nonliner_solve_picard' and the residual computation function had been
% changed to a parameter
minfo = electrical_network_init('R', 0.02, 'f0', 0.1, 'newton', false);
u0 = 0.01*rand(minfo.model_info.num_vars,1);
u2 = general_iterative_solver(minfo, p, 'u0', u0, 'verbose', true);

% compare both solutions (inline solver and function)
norm(u-u2)

%%
% now the same with Newton's method which should converge much faster
minfo = electrical_network_init('R', 0.02, 'f0', 0.1, 'newton', true);
u0 = 0.01*rand(minfo.model_info.num_vars,1);
u3 = general_iterative_solver(minfo, p, 'u0', u0, 'verbose', true);

norm(model_residual(minfo, u2, p))
norm(model_residual(minfo, u3, p))

