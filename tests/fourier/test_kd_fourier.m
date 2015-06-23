clf

a=2; b=3.5;
K=67;
%func = @(x)(exp(-abs(x/2)));
%func = @(x)(exp(-abs(x/1)));
func = @(x)(exp(-abs(x/10).^2));
[sigma_k, wp_k] = kd_fourier(func, [a, b], K);

xi = linspace(a,b,100);
[x1,x2]=meshgrid(xi);
C_ex=func(x1-x2);

r_k = sin_basis_eval(wp_k, xi);
u_k = binfun(@times, sigma_k, r_k);
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