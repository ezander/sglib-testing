function [u, a, pos]=diffusion_1d_complete_solve(a1, a2)

model=diffusion_1d_init();
[u,model]=model_solve(model, [a1;a2]);
a = a1 * model.a{1} + a2 * model.a{2};
pos = model.pos;
return
