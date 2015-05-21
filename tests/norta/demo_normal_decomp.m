% Set the number of Gaussians
m = 3;

% Set up a covariance matrix
A = rand(m, m);
C = A*A';
C = 0.5 * (C + C');
assert(norm(C-C')<1e-10) % make sure C is symmetric
assert(all(eig(C))>0) % make sure C is positive

% Do a Cholesky factorisation
L = chol(C, 'lower'); % need to use 'lower' (default is right upper matrix)
assert( norm(C - L*L' )<1e-10); % make sure L really is a factorisation of C

% Now generate N independent Gaussian samples
N = 10000;
xi = randn(m, N);

% The samples in g are now Gaussian but should have correction matrix C
g = L * xi;

% We can check that using the builtin cov function
C_g = cov(g');

% Now compare C and C_g (the differenve should be within the MC error),
% which is approximately norm(C)/sqrt(N)
strvarexpand('Exact cov matrix= $C$');
strvarexpand('Sample cov matrix=$C_g$');
strvarexpand('rel_err  = $norm(C_g-C)/norm(C)$');
strvarexpand('expected = $norm(C)/sqrt(N)$');




