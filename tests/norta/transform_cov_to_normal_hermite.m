function rho_N = transform_cov_to_normal_hermite(X1, X2, rho_X, M)
% TRANSFORM_COV_TO_NORMAL_HERMITE Compute covariance of Gaussians using Hermite expansion.

% Setup the polynomial (up to order M)
p = setup_polynomial(X1, X2, M, rho_X);

% Compute the roots of p
rs = roots(p);

% Only take roots between approx -1 and 1 and with imaginary part approx 0
delta_re = 1e-8;
delta_im = 1e-10;
ind=(real(rs)>=-1-delta_re) & ...
    (real(rs)<=1+delta_re) & ...
    (abs(imag(rs))<delta_im);

% Check that those roots exist and are unique
rho_N=unique(rs(ind));
if isempty(rho_N)
    error( 'the inverse of the transform polynomial could not be found.' );
elseif length(rho_N)>1
    error( 'the inverse of the transform polynomial is not unique.' );
end


function p_i = setup_polynomial(X1, X2, M, rho_X)
% SETUP_POLYNOMIAL Setup the polynomials for the Hermite expansion method

% Set the number of nodes for integration
n = min(M, 7);

% Compute the Hermite coefficients for both random vars X1 and X2
a_i = hermite_expand(X1, M, n);
b_i = hermite_expand(X2, M, n);

% Compute the expression given in the lecture (setting p_i(1) to zero means
% the same as subtracting the mean, because that's the coefficient of H_0)
p_i = a_i .* b_i .* factorial((0:M)');
p_i(1) = 0;

% Since we want to solve p(rho_N)=rho_X we have to subtract rho_X, so
% that we can solve p'(rho_N)=p(rho_N)-rho_X=0
p_i(1) = p_i(1) - rho_X;

% Matlab orders the polynomials cofficients from highest to lowest so we
% need to reorder them
p_i = p_i(end:-1:1)';
