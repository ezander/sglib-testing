function test2

N=20;
sig = 1;
b = 1;

[pos, els]=create_mesh_1d(0, 1, 100);

[ra_i_k, sigmaa_k] = kl_evp_analytic(pos, els, N, sig, b)
[rn_i_k, sigman_k] = kl_evp_numerically(pos, els, N, sig, b)

format short g
[sigmaa_k; sigman_k; sigmaa_k-sigman_k]'
1;

