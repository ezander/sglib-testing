clf

a=2; b=3.5;
%K=67;
%func = @(x)(exp(-abs(x/2)));
%func = @(x)(exp(-abs(x/1)));
l_c = 0.2;
%l_c = 1.5;
cov_func = @(x)(exp(-abs(x/l_c).^2));
pow_func = @(w)(l_c*sqrt(pi)*exp(-abs(pi*w*l_c).^2));

cov_func = @(x)(exp(-abs(x/l_c)));
pow_func = @(w)(2*l_c./(1+(2*pi*w*l_c).^2));



enlarge=true;
[sigma_k, wp_k] = kd_fourier(cov_func, [b-a], 'autoenlarge', enlarge);
[sigma_k, wp_k] = kd_fourier(pow_func, [b-a], 'is_spectral', true, 'autoenlarge', enlarge);

xi = linspace(a,b,100);
[x1,x2]=meshgrid(xi);
C_ex=cov_func(x1-x2);

u_k = trig_basis_eval(sigma_k, wp_k, xi);
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
