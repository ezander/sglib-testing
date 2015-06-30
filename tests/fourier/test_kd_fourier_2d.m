%%
clf

a=2; b=3.5;
l_c = 0.1;

nu = 0.5;
%nu = 2.0;
cov_func = @(x)(matern_covariance(nu, x,[],l_c));
pow_func = @(w, d)(matern_spectral_density(nu, w, l_c, d));

[pos,els]=load_pdetool_geom( 'lshape', 'numrefine', 3 );


enlarge=true;
options = {'autoenlarge', enlarge, 'max_funcs', 2000, 'ratio', 0.99};
%[sigma_k, wp_k] = kd_fourier(cov_func, pos, options);
[sigma_k, wp_k, u_k] = kd_fourier(pow_func, pos, 'is_spectral', true, options{:});

clf
for i=1:min(size(u_k,1),20)
    [i, sigma_k(i), wp_k(i,[1,3])]
    plot_field(pos, els, u_k(i,:)', 'show_mesh', false)
    drawnow
    pause(0.15);
end


for alpha=[linspace(0,pi/2), pi/4]
    x = [-1;-1] + 1.7 * [cos(alpha); sin(alpha)];
    u_x = trig_basis_eval(sigma_k, wp_k, x);
    plot_field(pos, els, u_k' * u_x, 'show_mesh', false)
    %view(3)
    drawnow
    pause(0.03);
end
