function [un, step_info, state] = diffusion_1d_step(state, u, a)

[un, step_info, state] = diffusion_1d_solve(state, a);
un = u + state.step_relax * (un - u);
