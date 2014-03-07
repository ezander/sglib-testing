subplot(1,2,1); 
[u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
hold on;

xi_plot = gpcgerm_sample(V_u, 5, 'mode', 'qmc');
a_plot = gpc_evaluate(a_i_alpha, V_a, xi_plot);
for i=1:5
    plot(pos, gpc_evaluate(u_i_alpha, V_u, xi_plot(:,i)), 'k:' );
    [u_ex, model] = model_solve(model, a_plot(:,i));
    plot(pos, u_ex, 'r:');
end
legend('mean-std', 'mean', 'mean+std', 'surrogate samples', 'exact samples'); 
title('resp. surf. proj.'); ylim([0,3.5]); grid on; hold off;
enhance_plot;


plot_response_surface(u_i_alpha([10,30,60,90],:), V_u); 
title('response surfaces at 0.1, 0.3, 0.6 and 0.9'); zlim([0, 5]);
enhance_plot;
