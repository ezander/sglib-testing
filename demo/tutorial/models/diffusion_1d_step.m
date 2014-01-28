function [du, state] = diffusion_1d_step(state, u, a)

[u_sol, solve_info, state] = diffusion_1d_solve(state, a);
du = 0.20 * (u_sol - u);
