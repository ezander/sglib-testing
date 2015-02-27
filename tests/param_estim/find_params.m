function [ampl_est, alpha_est, omega_est, phi_est] = find_params(x, y, do_opt, n)

xn = x(n);
x = x(n:end) - xn;
y = y(n:end);

[zero_ind, zero_dir] = zero_find(y);
peak_ind = peak_find(y);
zero_loc = (zero_ind-1) * (x(2)-x(1));

i=0:(length(zero_ind)-1);
zero_fit=polyfit(i, zero_loc, 1);

dt = zero_fit(1);
omega_est = pi / dt;

phi_est = -omega_est * zero_fit(2);
if zero_dir(1)==-1; phi_est = phi_est + pi; end

j=0:(length(peak_ind)-1);
yp = abs(y(round(peak_ind)));
peak_fit=polyfit(j*dt, -log(yp), 1);
alpha_est = peak_fit(1);

y_est = exp(-alpha_est*x).*sin(omega_est*x+phi_est);
a = polyfit(y_est, y, 1);
ampl_est = a(1);

if do_opt
    func = @(p)(optimfunc(x, y, p(1), p(2), p(3), p(4)));
    p0=[ampl_est; alpha_est; omega_est; phi_est];
    H0=eye(4);
    newton_opts.abstol = 1e-7;
    newton_opts.verbosity = 1;
    [p,flag,iter] = minfind_quasi_newton(func, p0, H0, newton_opts)
    
    ampl_est = p(1);
    alpha_est = p(2);
    omega_est = p(3);
    phi_est = p(4);
end

phi_est = phi_est - xn*omega_est;
ampl_est = ampl_est * exp(alpha_est * xn);
