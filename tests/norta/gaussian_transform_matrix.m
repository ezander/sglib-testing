function L=gaussian_transform_matrix(rho_N)
% GAUSSIAN_TRANSFORM_MATRIX Set up a transform matrix to compute correlated Gaussians.

C = [1, rho_N; rho_N, 1];
L = chol(C, 'lower');
