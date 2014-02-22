function [un, model] = diffusion_1d_step(model, u, a)

[un, solve_info, model] = diffusion_1d_solve(model, a); %#ok<ASGLU>
un = u + model.step_relax * (un - u);
