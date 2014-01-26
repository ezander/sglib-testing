clc
clear


sys = 'h'
m = 3

V_p = {sys, multiindex(m, 1)};

h1_mean = 200;
h1_std  = 0.6;

b1_mean = 3;
b1_std  = 0.2;

b2_mean = 100;
b2_std  = 5;

p_i_alpha = [
 h1_mean, h1_std, 0, 0;
 b1_mean,      0, b1_std, 0;
 b2_mean,      0, 0, b2_std;
 ]

p_i_alpha = [[h1_mean; b1_mean; b2_mean],  diag([h1_std, b1_std, b2_std])];


%% Monte Carlo
N = 10000;
theta_i = gpcgerm_sample(V_p, N);
p_i = gpc_evaluate(p_i_alpha, V_p, theta_i);
x_i = beam(p_i);

plot(p_i(1,:), p_i(2,:), '.')
plot3(p_i(1,:), p_i(2,:), p_i(3,:), '.')
axis equal
grid on

subplot(2,1,1)
empirical_density(x_i(1,:));

mean(x_i,2)


nio = x_i>5e4;
x_nio = x_i(nio);

subplot(2,1,2)
empirical_density(x_nio(1,:));
kde(x_i(1,:));

%%

x_mean = gpc_integrate(@beam, V_p, 2, 'grid', 'smolyak')

%%

V_p = {'Plh', multiindex(m, 1)};
[p, w] = gpc_integrate([], V_p, [4, 7, 20],  'grid', 'full_tensor');
plot3(p(1,:), p(2,:), p(3,:), '.')

%%
[p, w] = gpc_integrate([], V_p, 4, 'gpc_coeffs', p_i_alpha, 'grid', 'smolyak');
p
plot(p(1,:), p(2,:), '.')

N=length(w);

wmax = 5e4;
nio=0;
for i=1:N
    delta = beam(p(:,i));
    nio = nio + (delta>wmax) * w(i);
end
nio

