%%
clear
clc
clf
multiplot_init(1, [], 'ordering', 'row')

%% Set up the model
q_s = 1.1;
q_m = 10*[0.1; 0.3];

% The solution opertor q_s * (q - q_m)^2
S = @(q)( q_s * sum((binfun(@minus, q, q_m)).^2, 1) );
% The measurement operator
M = @(u, q)( u );
% The combined measurement and solution operator
G = @(q)( M(S(q), q) );

% The true parameter
q_true = 3*[0.4; 0.2];
q_true = q_m;
q_true = [0.0; 0.0];
% The true measurement value at q_true
y_true = G(q_true);
% The measured value (currently assumed to be the true value)
y_m = y_true;
var_factor = 1;
sigma_eps = var_factor * 0.1;


%% Create a surrogate model for the parameters
V_Q0 = gpcbasis_create('H', 'm', size(q_true, 1), 'p', 2);
Q0_i_alpha = gpc_projection(@(xi)(xi), V_Q0, 2);
Y_func = @(xi)G(gpc_evaluate(Q0_i_alpha, V_Q0, xi));

Y_i_alpha = gpc_projection(Y_func, V_Q0);

%%
%for p_phi=1:6
p_phi=3;

E=generate_stdrn_simparamset(sigma_eps);
E_func = funcreate(@germ2params, E, @funarg);
V_E = E.get_germ();

%p_phi=2;
p_pn=2*p_phi+1;
p_int_mmse=20;
p_int_proj=20;
[Q1_i_beta, V_Q1] = mmse_update_gpc(Q0_i_alpha, Y_func, V_Q0, y_m, E_func, V_E, p_phi, p_int_mmse, p_pn, p_int_proj, 'int_grid', 'full_tensor');

%% Plot parameter prior and posterior pdf
multiplot
Ns=30000;
Q0_samples = gpc_sample(Q0_i_alpha, V_Q0, Ns);
Q1_samples = gpc_sample(Q1_i_beta, V_Q1, Ns);

%%
mcmc_func = @(xi)( E.pdf(y_m - G(gpc_evaluate(Q0_i_alpha, V_Q0, xi))) .* gpcgerm_pdf(V_Q0, xi) );
[pos,els]=create_mesh_2d_rect(8);
pos=4*pos-2;
z=funcall(mcmc_func, pos)';
z=z/max(z(:));

%plot_field(pos, els, 10*z-0.01, 'show_mesh', false)


%%
N=1000;
N=Ns;
prop_dist = GPCGermDistribution(V_Q0);
xi0 = gpcgerm_moments(V_Q0);
QMCMC_samples = mh_sample_parallel(N, xi0, mcmc_func, prop_dist);


%%
hold off;
dY_i_alpha = Y_i_alpha;
dY_i_alpha(:,1) = dY_i_alpha(:,1) - y_m;

plot_response_surface(dY_i_alpha, V_Q0, 'pdf_plane', 0, 'delta', 0.001, 'N', 40, 'surf_color', 0.25);
hold all;

zzz = zeros(Ns,1);
plot3(Q1_samples(1,:), Q1_samples(2,:), zzz+0.02, 'g.', 'MarkerSize', 3)
plot3(QMCMC_samples(1,:), QMCMC_samples(2,:), zzz+0.04, 'r.', 'MarkerSize', 5)
zlim([-1,1])

h=plot_field(pos, els, 0.5*z-0.001, 'show_mesh', false);

axes1 = gca;
%view(axes1,[50 22]);
view(axes1,[63.2 26.8]);
grid(axes1,'on');
axis(axes1,'square');

save_png( gcf, strvarexpand('compare_mmse_mcmc_$var_factor$'), 'figdir', '.')
