%%
% initialise the model of the electrical network
clc
model = electrical_network_init('R', 0.02, 'f0', 0.1);
% display the model
model

% get the number of parameters from the model and set to some random value
p = rand(model.model_info.num_params, 1);
% get a starting vector from the model 
u = funcall(model.model_info.sol_init_func, p);
% compute the residual for the starting vector and the given model
res = electrical_network_residual(model, u, p);
norm(res)

%%
% solve the system with a nonlinear solver based on Picard iterations
model = electrical_network_init('R', 0.02, 'f0', 0.1, 'newton', false);
u = funcall(model.model_info.sol_init_func, p);
maxiter = 100;
abstol = 1e-5;

for iter = 1:maxiter
    res = model_residual(model, u, p);
    fprintf( 'iter %d, norm=%g\n', iter, norm(res));
    if norm(res)<abstol
        break;
    elseif iter==maxiter
        error('solve:no_conv', ...
            'Could not reach convergence (abstol=%g) in %d iterations', ...
            abstol, maxiter);
    end
    
    u = model_step(model, u, p);
end

%%
% the same as before, but solver has been put into the function
% 'nonliner_solve_picard' and the residual computation function had been
% changed to a parameter
model = electrical_network_init('R', 0.02, 'f0', 0.1, 'newton', false);
u0 = 0.01*rand(model.model_info.num_vars,1);
u2 = general_iterative_solver(model, p, 'u0', u0, 'verbosity', 1);

% compare both solutions (inline solver and function)
norm(u-u2)

%%
% now the same with Newton's method which should converge much faster
model = electrical_network_init('R', 0.02, 'f0', 0.1, 'newton', true);
u0 = 0.01*rand(model.model_info.num_vars,1);
u3 = general_iterative_solver(model, p, 'u0', u0, 'verbosity', 1);

norm(model_residual(model, u2, p))
norm(model_residual(model, u3, p))

