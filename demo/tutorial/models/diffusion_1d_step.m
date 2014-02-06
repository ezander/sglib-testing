function [un, state] = diffusion_1d_step(state, u, a)

[un, solve_info, state] = diffusion_1d_solve(state, a); %#ok<ASGLU>
un = u + state.step_relax * (un - u);
