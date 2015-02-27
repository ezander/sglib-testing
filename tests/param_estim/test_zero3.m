function test_zero3

clc
clf

alpha1 = 0.2;
omega1 = 1.4;
phi1 = 0.123;

alpha2 = 1.6;
omega2 = 3.4;
phi2 = 0.45;

x = linspace(0, 10, 100);
y1 = exp(-alpha1*x).*sin(omega1*x + phi1);
y2 = -0.6 * exp(-alpha2*x).*sin(omega2*x + phi2);
y = y1 + y2;

n = 30;
[ampl1_est, alpha1_est, omega1_est, phi1_est] = find_params(x, y, true, n);
y1_est = ampl1_est * exp(-alpha1_est*x).*sin(omega1_est*x+phi1_est);

[ampl2_est, alpha2_est, omega2_est, phi2_est] = find_params(x, y-y1_est, true, 1);
y2_est = ampl2_est * exp(-alpha2_est*x).*sin(omega2_est*x+phi2_est);

y_est = y1_est + y2_est;





plot(x, y)
hold all
plot(x, y_est+0.001);
plot(x, y-y_est-0.001);
legend('y', 'y_est', 'dy');
hold off;


norm(y-y_est)

