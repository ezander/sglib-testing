%%
clf

a=2; b=3.5;
l_c = 0.2;

nu = 1;
cov_func = @(x)(matern_covariance(nu, x,[],l_c));
pow_func = @(w, d)(matern_spectral_density(nu, w, l_c, d));

[pos,els]=load_pdetool_geom( 'lshape', 'numrefine', 1 );


enlarge=true;
options = {'autoenlarge', enlarge, 'max_funcs', 300, 'ratio', 0.99};
%[sigma_k, wp_k] = kd_fourier(cov_func, pos, options);
[sigma_k, wp_k, u_k] = kd_fourier(pow_func, pos, 'is_spectral', true, options{:});

clf
for i=1:min(size(u_k,1),10)
    [i, wp_k(i,[1,3])]
    plot_field(pos, els, u_k(i,:)')
    drawnow
    pause(0.1);
end
plot_field(pos, els, sum(u_k,1)')


%% 

%
if false
xi = linspace(a,b,100);
[x1,x2]=meshgrid(xi);
C_ex=reshape(cov_func(x1(:)'-x2(:)'), size(x1));

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
end
