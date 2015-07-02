clear
l_c = 0.02;

cov_func = @(x)(exponential_covariance(x,[],l_c));
pow_func = @(w,d)(exponential_spectral_density(w,l_c,d));

cov_func = @(x)(gaussian_covariance(x,[],l_c));
pow_func = @(w,d)(gaussian_spectral_density(w,l_c,d));

nu = 1.5;
cov_func = @(x)(matern_covariance(nu, x,[],l_c));
pow_func = @(w,d)(matern_spectral_density(nu, w,l_c,d));



N = 301;
L = 3;
W = pi / 2 / L * N;

[X, Y] = meshgrid(linspace(-L,L,N)); x = [X(:),Y(:)]';
[WX, WY] = meshgrid(linspace(-W,W,N)); w = [WX(:),WY(:)]';
sz = size(X);

C = reshape(cov_func(x), sz);
S = reshape(pow_func(w/(2*pi),2)/(2*L)^2, sz);

multiplot_init(2,2,'ordering','row');
surf_options={'FaceColor','interp','EdgeColor','none','FaceLighting','phong'};

multiplot;
surf(X,Y,C,surf_options{:})
multiplot;
C2 = ifft2(S);
C2 = abs(C2)*N^2;
C2 = fftshift(C2);
surf(X,Y,C2,surf_options{:})

multiplot;
surf(WX,WY,S,surf_options{:})
multiplot;
S2 = fft2(C);
S2 = abs(S2)/N^2;
S2 = fftshift(S2);
surf(WX,WY,S2,surf_options{:})

multiplot_adjust_range('separate', 'rows');

