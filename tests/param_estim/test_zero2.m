function test_zero2

clc
clf

alpha = 0.2;
omega = 1.4;
phi = 0.123;

alpha2 = 1.6;
omega2 = 3.4;
phi2 = 0.45;

x = linspace(0, 10, 100);
y = exp(-alpha*x).*sin(omega*x + phi);
y2 = exp(-alpha2*x).*sin(omega2*x + phi2);
y = y + y2;

n = 36;
[ampl_est, alpha_est, omega_est, phi_est] = find_params(x, y, true, n);
y_est = ampl_est * exp(-alpha_est*x).*sin(omega_est*x+phi_est);

plot(x, y)
hold all
plot(x, y_est+0.001);
plot(x, y-y_est-0.001);
legend('y', 'y_est', 'dy');
hold off;


norm(y-y_est)

