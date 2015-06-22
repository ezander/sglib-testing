function [y_k_i] = sin_basis_eval(sin_wp_k, x_i)
w_k=sin_wp_k(:,1); 
p_k=sin_wp_k(:,2); 
y_k_i = sin(binfun(@plus, w_k*x_i, p_k));

