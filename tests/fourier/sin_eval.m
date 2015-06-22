function [y_j_i] = sin_eval(A_j_k, sin_wp_k, x_i)
y_k_i = sin_basis_eval(sin_wp_k, x_i);
y_j_i = A_j_k * y_k_i;


