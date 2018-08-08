sys = RealSystem()

clf
hold on
x_true_i = [];
y_true_i = [];
u_true_i = [];
v_true_i = [];
x_est_i = [];
y_est_i = [];
u_est_i = [];
v_est_i = [];


x0 = [0.1; 0.2; 1.2; 0];
P0 = eye(4) * 0.2;

xp = x0;
Pp = P0;

F = sys.get_transition_matrix();
H = sys.get_observation_matrix();
Q = sys.get_process_noise_matrix();
R = sys.get_observation_noise_matrix();

for i=1:1000
    sys.do_step()
    x_true_i(end+1) = sys.x;
    y_true_i(end+1) = sys.y;
    u_true_i(end+1) = sys.u;
    v_true_i(end+1) = sys.v;
    
    
    xnp = F * xp;
    Pnp = F * Pp * F' + Q;
    
    y = sys.get_observation();
    S = H * Pnp * H' + R;
    K = Pnp * H' * inv(S);
    xn = xnp + K * (y - H*xnp);
    I = eye(4);
    Pn = (I - K*H) * Pnp * (I - K*H)' + K * R * K';
    
    x_est_i(end+1) = xn(1);
    y_est_i(end+1) = xn(2);
    u_est_i(end+1) = xn(3);
    v_est_i(end+1) = xn(4);
    subplot(2,1,1)
    plot(x_true_i, y_true_i, x_est_i, y_est_i,'-k')
    subplot(2,1,2)
    plot(u_true_i, v_true_i, u_est_i, v_est_i,'-k')
    drawnow
    
    xp = xn;
    Pp = Pn;
end
