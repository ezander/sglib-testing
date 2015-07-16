%% Show covariance
clear;
clf;

[cov_func, x, els] = setup_model();

[X,Y] = meshgrid(x);

C=covariance_matrix(x, cov_func);
surf(X,Y,full(C))
