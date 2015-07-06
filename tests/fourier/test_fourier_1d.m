clear
l_c = 0.2;

cov_func = @(x)(gaussian_covariance(x,[],l_c));
pow_func = @(w,d)(gaussian_spectral_density(w,l_c,d));

nu = 1.5;
cov_func = @(x)(matern_covariance(nu, x,[],l_c));
pow_func = @(w,d)(matern_spectral_density(nu, w,l_c,d));

cov_func = @(x)(exponential_covariance(x,[],l_c));
pow_func = @(w,d)(exponential_spectral_density(w,l_c,d));


N = 301;
L = 3;
W = pi / 2 / L * N;

x = linspace(-L,L,N);
w = linspace(-W,W,N);
sz = size(x);

C = reshape(cov_func(x), sz);
S = reshape(pow_func(w/(2*pi),1)/(2*L), sz);

multiplot_init(2,1);
plot_options={};

multiplot;
plot(x,C); 
hold all;
C2 = ifft2(S);
C2 = abs(C2);
C2 = fftshift(C2);
plot(x,C2)

multiplot;
plot(w,S+0.001); hold all;
S2 = fft2(C);
S2 = abs(S2)/N;
S2 = fftshift(S2);
plot(w,S2+0.002); hold all;

[A_k, TP, x_i] = fourier_series_expand(cov_func, -L, L, 300, 'symmetry', 'even');
A_k(2:end) = A_k(2:end)/2;
plot(2*pi*TP{1}, A_k+0.003)
