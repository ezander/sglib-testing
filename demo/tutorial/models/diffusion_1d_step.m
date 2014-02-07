function [un, minfo] = diffusion_1d_step(minfo, u, a)

[un, solve_info, minfo] = diffusion_1d_solve(minfo, a); %#ok<ASGLU>
un = u + minfo.step_relax * (un - u);
