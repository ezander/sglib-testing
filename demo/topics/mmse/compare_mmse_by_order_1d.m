%%
clear
clc
multiplot_init(6, [], 'ordering', 'row', 'separate_figs', false)

%% Set up the model
q_s = 0.03;
q_m = -3;

% The solution opertor
S = @(q)( q_s * (q - q_m).^2 );
% The measurement operator
M = @(u, q)( u );
% The combined measurement and solution operator
G = @(q)( M(S(q), q) );

% The true parameter
q_true = 3;
% The true measurement value at q_true
y_true = G(q_true);
% The measured value (currently assumed to be the true value)
y_m = y_true;
sigma_eps = 2*0.7*q_s;


%% Create a surrogate model for the parameters
V_Q0 = gpcbasis_create('H', 'p', 1);
Q0_i_alpha = gpc_projection(@(xi)(2*xi), V_Q0, 2);

%%
for p_phi=1:6
    Y_func = @(xi)G(gpc_evaluate(Q0_i_alpha, V_Q0, xi));
    
    E=generate_stdrn_simparamset(sigma_eps);
    E_func = funcreate(@germ2params, E, @funarg);
    V_E = E.get_germ();
    
    %p_phi=2;
    p_pn=2*p_phi;
    p_int_mmse=20;
    p_int_proj=20;
    [Q1_i_beta, V_Q1] = mmse_update_gpc(Q0_i_alpha, Y_func, V_Q0, y_m, E_func, V_E, p_phi, p_int_mmse, p_pn, p_int_proj, 'int_grid', 'full_tensor');
    
    %% Plot parameter prior and posterior pdf
    multiplot
    Ns=100000;
    Q0_samples = gpc_sample(Q0_i_alpha, V_Q0, Ns);
    Q1_samples = gpc_sample(Q1_i_beta, V_Q1, Ns);
    plot_density(Q0_samples)
    hold all;
    plot_density(Q1_samples)
    Q_range = xlim;
    Q_range = linspace(Q_range(1), Q_range(2));
    plot(Q_range, G(Q_range))
    plot(Q_range, repmat(y_m + sigma_eps*[-1;0;1], 1, length(Q_range)))
    hold off;
    line([q_true, q_true], ylim, [-0.1,-0.1], 'LineWidth',2, 'Color',[.8 .8 .8])
    axis tight
end

multiplot_adjust_range('axes', 'x', 'range', [-5,5])
multiplot_adjust_range('axes', 'y', 'range', [0,2])
