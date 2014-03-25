function demo_mmse_for_ar

%%

m0=2; m=2; n=2; p_phi=1; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'exact2_p1');

m0=2; m=2; n=2; p_phi=2; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'exact2_p2');

m0=2; m=2; n=2; p_phi=3; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'exact2_p3');


m0=2; m=1; n=4; p_phi=1; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'one_measurement_p1');

m0=2; m=1; n=4; p_phi=3; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'one_measurement_p2');

m0=5; m=3; n=3; p_phi=3; p_int=7; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'hidden5');


% m0=3; m=2; n=2; p_phi=3; p_int=9; p_x=3; p_y=2;
% demo(m0, m, n, p_phi, p_int, p_x, p_y, 'hidden_p3');




function demo(m0, m, n, p_phi, p_int, p_x, p_y, str)
rand_seed(935590374)
V_y=gpcbasis_create('H', 'm', m0, 'p', p_y);
y_j_beta = gpc_rand_coeffs(V_y, m);

V_x=gpcbasis_create('H', 'm', m0, 'p', p_x);
x_i_alpha = gpc_rand_coeffs(V_x, n);

[phi_j_gamma, V_phi]=mmse_estimate(x_i_alpha, V_x, y_j_beta, V_y, p_phi, p_int, 'cond_warning', 1e10);

x_func = funcreate(@gpc_evaluate, x_i_alpha, V_x, @funarg);
y_func = funcreate(@gpc_evaluate, y_j_beta, V_y, @funarg);

[phi2_j_gamma, V_phi2]=mmse_estimate_generic(x_func, y_func, V_y, p_phi, p_int, 'cond_warning', 1e10);

assert_equals(phi_j_gamma, phi2_j_gamma);
assert_equals(V_phi, V_phi2);




%%
N=100;
xi = gpcgerm_sample(V_y,N);
y = gpc_evaluate(y_j_beta, V_y, xi);
x = gpc_evaluate(x_i_alpha, V_x, xi);
x_approx = gpc_evaluate(phi_j_gamma, V_phi, y);

ms=10;
plot(x(1,:), x(2,:), 'ko', 'MarkerSize', ms); hold all;
plot(x_approx(1,:), x_approx(2,:), 'rx', 'MarkerSize', ms); hold off;
axis tight;
%return

if exist('str', 'var')
    disp(['saving: ', str, '.png ...']);
    %print('-dpng', [str, '.png']);
    save_figure(gcf, str, 'figdir', jb_figdir);
end
