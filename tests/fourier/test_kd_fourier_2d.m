%%
clf

l_c = 0.2;

nu = 0.5;
nu = 2.0;
cov_func = @(x)(matern_covariance(nu, x,[],l_c));
pow_func = @(w, d)(matern_spectral_density(nu, w, l_c, d));

cov_func = @(x)(gaussian_covariance(x,[],l_c));
pow_func = @(w,d)(gaussian_spectral_density(w,l_c,d));

cov_func = @(x)(exponential_covariance(x,[],l_c));
pow_func = @(w,d)(exponential_spectral_density(w,l_c,d));

[pos,els]=load_pdetool_geom( 'lshape', 'numrefine', 3 );


enlarge=true;
options = {'autoenlarge', enlarge, 'max_funcs', 2700, 'ratio', 0.99};
%[sigma_k, wp_k] = kd_fourier(cov_func, pos, options);
[sigma_k, wp_k, u_k] = kd_fourier(pow_func, pos, 'is_spectral', true, options{:});

clf
multiplot_init(2,3);

multiplot
plot(wp_k(:,1), wp_k(:,3), '.'); axis equal;

multiplot;
x = [0.0; 0.0];
w = linspace(-1,1, 1001);
u_x = trig_basis_eval(sigma_k, wp_k, x);

plot(w, cov_func([1;0]*w)); hold all;
legend_add('exact');
for alpha = linspace(0, pi/2, 5)
    y = [cos(alpha); sin(alpha)] * w;
    u_y = trig_basis_eval(sigma_k, wp_k, y);
    plot(w, u_y' * u_x + 0*alpha*0.03);
    legend_add(alpha);
end
hold off;


multiplot;
for i=1:min(size(u_k,1),0*20)
    [i, sigma_k(i), wp_k(i,[1,3])]
    plot_field(pos, els, u_k(i,:)', 'show_mesh', false)
    drawnow
    pause(0.15);
end


h = multiplot;
for alpha=[linspace(0,pi/2,40), pi/4]
    x = [-1;-1] + 1.7 * [cos(alpha); sin(alpha)];
    u_x = trig_basis_eval(sigma_k, wp_k, x);
    multiplot(h)
    plot_field(pos, els, u_k' * u_x, 'show_mesh', false)
    view(3)
    multiplot
    plot_field(pos, els, u_k' * u_x, 'show_mesh', false)
    drawnow
    pause(0.02);
end
