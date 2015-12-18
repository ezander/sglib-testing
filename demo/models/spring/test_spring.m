%% Generate the model;
model = SpringModel('T', 20);

%% Generate the parameter set for the prior model
Q = SimParamSet();
Q.add('u0', 1)
Q.add('v0', 1)
Q.add('m', 1)
Q.add('d', 0.03)
Q.add('k', 1)

Q.set_dist('m', SemiCircleDistribution().fix_bounds(0.5, 2.5));
Q.set_dist('k', SemiCircleDistribution().fix_bounds(0.5, 2.5));

%% Generate some artificial truth we want to identify which is (for the
% moment) close to the mean of the prior model
g_func = @(q)(model.compute_measurements(model.compute_response(q)));
q_true = mean(Q.sample(4),2);
y_true = g_func(q_true);
y_m = y_true + normal_sample(1, 0, 0.01);

% Generate measurement function
Y_func = @(xi)(compute_measurement(xi, g_func, Q));

%% Generate surrogate model and replace original model with it
surr_model = generate_surrogate_model(model, Q, 5, 'projection', {5, 'grid', 'full_tensor'});
surr_model.coeffs

u_i_alpha = surr_model.coeffs;
V_u = surr_model.basis;
plot_response_surface(u_i_alpha(1,:), V_u);

% model = surr_model

%% Generate the error model
sigma_eps = 0.2;
E = SimParamSet();
for i=1:length(y_m)
    param_name=strvarexpand('E_y$i$');
    E.add(param_name, NormalDistribution(0,sigma_eps));
end
% Get error function and gpc germ
E_func = @(xi)(E.germ2params(xi));
V_e = E.get_germ();



%% Generate original parameter GPC and update with the MMSE
[Q_i_alpha, V_q] = Q.gpc_expand();
p_phi=1;
p_pn=2;
p_int_mmse=3;
p_int_proj=4;
[Qn_i_beta, V_qn] = mmse_update_gpc(Q_i_alpha, Y_func, V_q, y_m, E_func, V_e, p_phi, p_int_mmse, p_pn, p_int_proj);


%% Compute some stuff with updated parameters
[q_mean, q_var] = gpc_moments(Q_i_alpha, V_q);
[qn_mean, qn_var] = gpc_moments(Qn_i_beta, V_qn);
gpc_covariance(Q_i_alpha, V_q)
chopabs(gpc_covariance(Qn_i_beta, V_qn, [], 'corrcoeffs', true))
chopabs([q_mean-q_true, sqrt(q_var), qn_mean-q_true, sqrt(qn_var)])
chopabs([g_func(q_mean)-y_m, g_func(qn_mean)-y_m, y_true-y_m ])

