%% Compute KL Eigenfunctions
clear;
clf;

[cov_func, x, els, l_c] = setup_model();

% Compute covariance matrix
C=covariance_matrix(x, cov_func);
C = full(C);

% Compute Gramian matrix (aka mass matrix)
G = mass_matrix(x, els);
G = full(G);

% Solve the KL eigenvalue problem
[V,D] = eig(G*C*G, G);

% Take only KL eigenvalues larger 10^-14
lambda = diag(D);
ind = lambda>1e-14;
lambda = lambda(ind);

% Select and normalise the KL eigenfunctions
R = V(:,ind);
R_norm = sqrt(diag(R'*G*R)');
R = binfun(@times, R, 1./R_norm)

% Compute the sigma's
sigma = sqrt(lambda);


multiplot_init(2,2);

% Plot the numerically computed eigenvalues
multiplot;
semilogy(sigma, 'x')

% Plot the numerically computed eigenfunctions
multiplot;
plot(x, R(:,1:5))


% Compute the KL eigenvalues and functions analytically
N = length(sigma);
[R2, sigma2] = kl_solve_1d_exp(x, 1, l_c, N)
% Uncomment the next line to have the comparison KL computed numerically by
% the standard sglib function
% [R2, sigma2] = kl_solve_evp(C, G, N)



% Plot the analytical eigenvalues
multiplot;
semilogy(sigma2, 'x')

% Plot the analytical eigenvalues
multiplot;
R2=flip_to_align(R2, R);
plot(x, R2(:,1:5))
