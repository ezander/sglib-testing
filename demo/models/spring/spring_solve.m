function [x,v] = spring_solve(u0, v0, m, d, k, T, varargin)

options = varargin2options(varargin);
[solver, options] = get_option(options, 'solver', 'direct');
check_unsupported_options(options, mfilename);

switch solver
    case 'direct'
        [x,v]=undamped_spring_solve_direct(u0, v0, m, k, d, T);
    case 'numerical'
        [x, v]=undamped_spring_solve_ode(u0, v0, m, k, d, T);
    otherwise
        error('foo:bar', 'Unknown solver: %s', solver);
end

function [xt, vt]=undamped_spring_solve_direct(x0, v0, m, k, d, T)
assert(k>0 && m>0 && 0*d>=0);
alpha = d/m;
D = k/m-(d/m)^2;
omega = sqrt(D);
xt = exp(-alpha*T) * (x0 * cos(omega*T) + (v0 + alpha*x0) / omega * sin(omega*T));
vt = exp(-alpha*T) * (v0 * cos(omega*T) - ((alpha*v0 + alpha^2*x0) / omega + x0 * omega) * sin(omega*T));


function [xt, vt, t, ut]=undamped_spring_solve_ode(u0, m, k, d, T)
[t,ut] = ode45(@undamped_spring_ode, [0, T], u0, [], m, k, d);
assert(t(end)==T)
xt = ut(end,1);
vt = ut(end,2);

function [dudt] = undamped_spring_ode(~, u, m, k, d)
x = u(1);
v = u(2);
a = -(2*d*v+k*x)/m;
dudt = [v; a];
