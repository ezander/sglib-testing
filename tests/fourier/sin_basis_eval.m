function [y_k_i] = sin_basis_eval(sin_wp_k, x_i)
y_k_i = ones(size(sin_wp_k, 1), size(x_i, 2));
for j=1:size(x_i,1)
    w_k=sin_wp_k(:,1);
    p_k=sin_wp_k(:,2);
    y_k_i = y_k_i .* sin(binfun(@plus, w_k*x_i, p_k));
end
