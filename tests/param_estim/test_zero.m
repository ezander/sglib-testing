function test_zero

clc
clf

alpha = 0.2;
omega = 1.4;
phi = 0.123;

alpha2 = 1.6;
omega2 = 1.4;
phi2 = 0.45;

x = linspace(0, 10, 100);
y = exp(-alpha*x).*sin(omega*x + phi);
y2 = exp(-alpha2*x).*sin(omega2*x + phi2);
y = y + y2;

[zero_ind, zero_dir] = zero_find(y);
[peak_ind, peak_minmax] = peak_find(y);
zero_loc = (zero_ind-1) * (x(2)-x(1));
peak_loc = (peak_ind-1) * (x(2)-x(1));

i=0:(length(zero_ind)-1);
zero_fit=polyfit(i, zero_loc, 1);

dt = zero_fit(1);
omega_est = pi / dt;
[omega_est, omega]

phi_est = -omega_est * zero_fit(2);
if zero_dir(1)==-1; phi_est = phi_est + pi; end
[phi_est, phi]

j=0:(length(peak_ind)-1);
yp = abs(y(round(peak_ind)));
peak_fit=polyfit(j*dt, -log(yp), 1);
alpha_est = peak_fit(1);
[alpha_est, alpha]



subplot(2,2,1)
plot(x,y);
hold on
plot( zero_loc, 0*zero_ind, 'ro');
plot( peak_loc, y(round(peak_ind)), 'go');
hold off

subplot(2,2,2)
plot(i, zero_loc, 'x');
hold on
plot(i, zero_fit(2) + zero_fit(1) * i, 'r');
hold off

subplot(2,2,3)
plot(j*dt, -log(yp), 'x');
hold on
plot(j*dt, peak_fit(2) + peak_fit(1) * j * dt, 'g');
hold off

%%
subplot(2,2,4)

y_est = exp(-alpha_est*x).*sin(omega_est*x);
plot(x, y-y_est);


func = @(p)(optimfunc(x, y, p(1), p(2), p(3)));
p0=[alpha_est; omega_est; phi_est];
H0=eye(3);
newton_opts.abstol = 1e-7;
newton_opts.verbosity = 1;
[p,flag,iter] = minfind_quasi_newton(func, p0, H0, newton_opts)

alpha_est2 = p(1);
omega_est2 = p(2);
phi_est2 = p(3);
y_est2 = exp(-alpha_est2*x).*sin(omega_est2*x+phi_est2);
hold on;
plot(x, y-y_est2);
hold off;
format long g
[norm(y-y_est), norm(y-y_est2)]






