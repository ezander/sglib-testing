function [r_i_k, sigma_k] = kl_evp_numerically(pos, els, N, sig, b)
G_N=mass_matrix(pos,els);
C=covariance_matrix(pos,...
    funcreate(@exponential_covariance, @funarg, @funarg, b, sig))
[r_i_k,sigma_k]=kl_solve_evp( C, G_N, N)

