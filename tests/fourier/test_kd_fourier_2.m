clf

a=2; b=3.5;
l_c = 0.2;

cov_func = @(x)(exponential_covariance(x,[],l_c));
pow_func = @(w,d)(exponential_spectral_density(w,l_c,d));

enlarge=true;
options = {'autoenlarge', enlarge, 'max_funcs', 300, 'ratio', 0.99};
[sigma_k, wp_k] = kd_fourier(cov_func, [b-a], options);
[sigma_k, wp_k] = kd_fourier(pow_func, [b-a], 'is_spectral', true, options{:});

xi = linspace(a,b,100);
[x1,x2]=meshgrid(xi);

u_k = trig_basis_eval(sigma_k, wp_k, xi);
C = u_k.'*u_k;

clf

plot(xi,C(10,:));
