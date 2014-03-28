function [xn_i_beta, V_xn]=mmse_update_model(x_func, y_func, V_xy, ym_beta, V_ym, p_phi, p_int_mmse, p_un, p_int_proj)

[phi_j_delta,V_phi]=mmse_estimate(x_func, y_func, V_xy, p_phi, p_int_mmse);

phi_func = funcreate(@gpc_evaluate, phi_j_delta, V_phi);
ym_func = funcreate(@gpc_evaluate, ym_beta, V_ym);
xn_func = funcompose(ym_func, phi_func);

V_xn = gpcbasis_create(V_ym, 'p', p_un);
xn_i_beta = gpc_projection(xn_func, V_xn, p_int_proj);


