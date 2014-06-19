%% Init stuff

init_func = @electrical_network_init;
solve_func = @electrical_network_solve;
step_func = @electrical_network_picard_iter_step;

init_func = @undamped_spring_init;
solve_func = @undamped_spring_solve;
step_func = @foobar; %@undamped_spring_solve;
polysys = 'u';
init_func = @()(undamped_spring_init('T', 20, 'd', 0.1));
polysys = 'p';

%% Monte Carlo

N = 100;
[u_mean, u_var] = compute_moments_mc(init_func, solve_func, polysys, N);
show_mean_var('Monte-Carlo', u_mean, u_var)

%% Quasi Monte Carlo

[u_mean, u_var] = compute_moments_mc(init_func, solve_func, polysys, N, 'mode', 'qmc');
show_mean_var('Quasi Monte-Carlo', u_mean, u_var)

%% Latin hypercube

[u_mean, u_var] = compute_moments_mc(init_func, solve_func, polysys, N, 'mode', 'lhs');
show_mean_var('Latin hypercube', u_mean, u_var)

%% Direct integration full tensor grid

p = 5;
[u_mean, u_var] = compute_moments_quad(init_func, solve_func, polysys, p, 'grid', 'full_tensor');
show_mean_var('Full tensor grid integration', u_mean, u_var);

%% Direct integration sparse grid

p = 5;
[u_mean, u_var] = compute_moments_quad(init_func, solve_func, polysys, p, 'grid', 'smolyak');
show_mean_var('Sparse grid (Smolyak) integration', u_mean, u_var);




%% Projection with full tensor grid
p_u = 10;
p_int = p_u+1;

V_u = gpcbasis_create(polysys, 'm', 2, 'p', p_u);
[u_i_alpha, x, w] = compute_response_surface_projection(init_func, solve_func, V_u, p_int, 'grid', 'full_tensor');

[u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
show_mean_var('Projection (L_2, response surface, tensor)', u_mean, u_var);


% Plot the response surface
hold off;
plot_response_surface(u_i_alpha(1,:), V_u, 'delta', 0.01, 'N', 80);

u=gpc_evaluate(u_i_alpha, V_u, x);
hold on; plot3(x(1,:), x(2,:), u(1,:), 'rx'); hold off;

%% Projection with sparse grid
mh = multiplot_init(3,3)
for p_u = 1:9
    multiplot;
    p_int = p_u+1;
    
    
    V_u = gpcbasis_create(polysys, 'm', 2, 'p', p_u);
    [u_i_alpha, x, w] = compute_response_surface_projection(init_func, solve_func, V_u, p_int, 'grid', 'full_tensor');
    
    [u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
    show_mean_var('Projection (L_2, response surface, sparse)', u_mean, u_var);
    
    
    % Plot the response surface
    hold off;
    plot_response_surface(u_i_alpha(1,:), V_u, 'delta', 0.01);
    
    u=gpc_evaluate(u_i_alpha, V_u, x);
    hold on; plot3(x(1,:), x(2,:), u(1,:), 'rx'); hold off;
    drawnow;
    u_proj_i_alpha = u_i_alpha;
end
same_scaling(mh, 'z', 'range', [-2, 2]);
same_scaling(mh, 'c');
%plot_response_surface(u_i_alpha(1,:), V_u, 'delta', 0.01, 'surf_color', 'pdf', 'pdf_plane', 'none');
%plot_response_surface(u_i_alpha(1,:), V_u, 'delta', 0.01, 'surf_color', 'pdf');hold on
%plot_response_surface(u_i_alpha(2,:), V_u, 'delta', 0.01, 'surf_color', 'pdf');
%plot_response_surface(u_i_alpha(3,:), V_u, 'delta', 0.01, 'surf_color', 'pdf');hold off;


%% Full tensor grid collocation (interpolation)

p_u = 10;

V_u = gpcbasis_create(polysys, 'm', 2, 'p', p_u, 'full_tensor', true);
[u_i_alpha, x] = compute_response_surface_tensor_interpolate(init_func, solve_func, V_u, p_u);

%ind=(multiindex_order(V_u{2})>=3);
%u_i_alpha(:,ind)=0;

[u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
show_mean_var('Interpolation, tensor (response surface)', u_mean, u_var);

% Plot the response surface
hold off;
plot_response_surface(u_i_alpha(1,:), V_u, 'delta', 0.01, 'N', 80);

u=gpc_evaluate(u_i_alpha, V_u, x);
hold on; plot3(x(1,:), x(2,:), u(1,:), 'rx'); hold off;

u_tensorcoll_i_alpha = u_i_alpha;

%% Sparse grid collocation (regression)

%% Non-intrusive Galerkin
% p_u = 3;
% state = funcall(init_func);
% 
% p_int = 4;
% V_u = gpcbasis_create(polysys, 'm', state.num_params, 'p', p_u, 'full_tensor', false);
% [u_i_alpha,x,w]=compute_response_surface_nonintrusive_galerkin(init_func, step_func, V_u, p_int, 'grid', 'full_tensor');
% 
% % Plot the response surface
% hold off;
% plot_response_surface(u_i_alpha(1,:), V_u, 'delta', 0.01);
% 
% u=gpc_evaluate(u_i_alpha, V_u, x);
% hold on; plot3(x(1,:), x(2,:), u(1,:), 'rx'); hold off;
% 
% u_galerkin_i_alpha = u_i_alpha;
% 
% %%
% i=1;
% u_plot=[u_galerkin_i_alpha(i, :); u_proj_i_alpha(i, :)];
% u_plot(1,1)=u_plot(1,1)+0.001;
% plot_response_surface(u_plot, V_u)
