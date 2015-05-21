function a_i = hermite_expand(X, M, n)
% HERMITE_EXPAND Compute Hermite expansion up to specified degree

% Note: This function already uses the GPC functions, which is a bit over
% the top, but was easier to implement that way.

% This defines a GPC basis in Hermite polynomials (H), in one variable
% ('m', 1) up to degree M ('p', M)
V=gpcbasis_create('H', 'm', 1, 'p', M);

% Compute integration points for this basis (i.e. the 1d Gauss-Hermite
% rule with n points; note: the [] could be a function to integrate, but we
% want to get back the points themselves, which is indicates by the [])
[xi,w]=gpc_integrate([], V, n);

% Now compute a_i = sum_j H_i(xi_j) F_X^-1(Phi(xi_j)) w_j / E(H_i^2)
% with H_i(xi_j) being the i-th row of gpcbasis_evaluate(V, xi), 
% F_X^-1(Phi(xi_j)) being gendist_stdnor(xi, X) and the rest should be
% evident (except probably ordering of indices in the involved arrays).
a_i = binfun(@times, gpcbasis_evaluate(V, xi), gendist_stdnor(xi, X)) ...
    * w ./ gpcbasis_norm(V, 'sqrt', false);
