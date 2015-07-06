clear
clf

a=2; b=3.5;
l_c = 3.1;

nu = 0.1;
cov_func = @(x)(matern_covariance(nu, x,[],l_c));
pow_func = @(w,d)(matern_spectral_density(nu, w,l_c,d));


cov_func = @(x)(exponential_covariance(x,[],l_c));
pow_func = @(w,d)(exponential_spectral_density(w,l_c,d));

cov_func = @(x)(gaussian_covariance(x,[],l_c));
pow_func = @(w,d)(gaussian_spectral_density(w,l_c,d));


enlarge=true;
options = {'autoenlarge', enlarge, 'max_funcs', 500, 'ratio', 0.99};
TB = kd_fourier(cov_func, [b-a], options);
TB = kd_fourier(pow_func, [b-a], 'is_spectral', true, options{:});

xi = linspace(a,b,100);
[x1,x2]=meshgrid(xi);
C_ex=reshape(cov_func(x1(:)'-x2(:)'), size(x1));

u_k = trig_basis_eval(TB, xi);
C = u_k.'*u_k;

multiplot_init(2,2);

multiplot
options={...
    'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong'};
surf(x1,x2,C_ex,options{:});

multiplot
surf(x1,x2,C,options{:});

multiplot
surf(x1,x2,C_ex); hold all;
surf(x1,x2,C); hold off;
view(135,0)
zlim([0,1.1]);

multiplot
M=size(u_k,1);
N=3;
r_i = randn(N,M)*u_k;
plot(xi, r_i);
grid on;

multiplot_adjust_range('axes', 'z', 'cols', [1])
