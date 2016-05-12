function [cov_func, x, els, l_c] = setup_model
% SETUP_MODEL Setup the default model for the KL demos

% Setup the domain
a = 2;
b = 4;
[x,els,bnd]=create_mesh_1d(a, b, 30);


% Setup the covariance function
l_c=1.3;
base_cov_func = @exponential_covariance;
%base_cov_func = @gaussian_covariance;
cov_func = funcreate(base_cov_func, funarg, funarg, l_c);
