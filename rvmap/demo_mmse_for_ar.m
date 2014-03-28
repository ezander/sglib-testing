function demo_mmse_for_ar
% DEMO_MMSE_FOR_AR Demonstration of the MMSE methods.
%
% See also MMSE_ESTIMATE, MMSE_ESTIMATE_GPC

%   Elmar Zander
%   Copyright 2014, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

%%
clf

m0=2; m=2; n=2; p_phi=1; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'exact2_p1');

m0=2; m=2; n=2; p_phi=2; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'exact2_p2');

m0=2; m=2; n=2; p_phi=3; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'exact2_p3');


m0=2; m=1; n=4; p_phi=1; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'one_measurement_p1');

m0=2; m=1; n=4; p_phi=3; p_int=9; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'one_measurement_p2');

m0=5; m=3; n=3; p_phi=3; p_int=7; p_x=3; p_y=2;
demo(m0, m, n, p_phi, p_int, p_x, p_y, 'hidden5');


% m0=3; m=2; n=2; p_phi=3; p_int=9; p_x=3; p_y=2;
% demo(m0, m, n, p_phi, p_int, p_x, p_y, 'hidden_p3');




function demo(m0, m, n, p_phi, p_int, p_x, p_y, str)
rand_seed(935590374)
V_y=gpcbasis_create('H', 'm', m0, 'p', p_y);
y_j_beta = gpc_rand_coeffs(V_y, m);

V_x=gpcbasis_create('H', 'm', m0, 'p', p_x);
x_i_alpha = gpc_rand_coeffs(V_x, n);

[phi_j_gamma, V_phi]=mmse_estimate_gpc(x_i_alpha, V_x, y_j_beta, V_y, p_phi, p_int, 'cond_warning', 1e10);

x_func = gpc_function(x_i_alpha, V_x);
y_func = gpc_function(y_j_beta, V_y);

[phi2_j_gamma, V_phi2]=mmse_estimate(x_func, y_func, V_y, p_phi, p_int, 'cond_warning', 1e10);

assert_equals(phi_j_gamma, phi2_j_gamma);
assert_equals(V_phi, V_phi2);




N=100;
xi = gpcgerm_sample(V_y,N);
y = gpc_evaluate(y_j_beta, V_y, xi);
x = gpc_evaluate(x_i_alpha, V_x, xi);
x_approx = gpc_evaluate(phi_j_gamma, V_phi, y);

ms=10;
plot(x(1,:), x(2,:), 'ko', 'MarkerSize', ms); hold all;
plot(x_approx(1,:), x_approx(2,:), 'rx', 'MarkerSize', ms); hold off;
axis tight;

if exist('str', 'var')
    disp(['saving: ', str, '.png ...']);
    save_figure(gcf, str, 'figdir', jb_figdir);
end
