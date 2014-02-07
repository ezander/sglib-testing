function [u, a, pos]=diffusion_1d_complete_solve(a1, a2)

minfo=diffusion_1d_init();
[u,minfo]=model_solve(minfo, [a1;a2]);
a = a1 * minfo.a{1} + a2 * minfo.a{2};
pos = minfo.pos;
return
