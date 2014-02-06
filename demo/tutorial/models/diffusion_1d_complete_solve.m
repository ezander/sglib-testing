function [u, a, pos]=diffusion_1d_complete_solve(a1, a2)

state=diffusion_1d_init();
[u,state]=model_solve(state, [a1;a2]);
a = a1 * state.a{1} + a2 * state.a{2};
pos = state.pos;
return
