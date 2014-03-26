function demo_mmse_for_spde(varargin)
% DEMO_MMSE_FOR_SPDE Short description of demo_mmse_for_spde.
%   DEMO_MMSE_FOR_SPDE Long description of demo_mmse_for_spde.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example demo_mmse_for_spde">run</a>)
%
% See also

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

p_phi = 3;
p_int_mmse = 3;
p_un = 3;
p_int_proj = 6;



[u_i_alpha, V_u, pos] = load_spde_solution('u');
[k_i_alpha, V_k] = load_spde_solution('k');

u_func = funcreate(@gpc_evaluate, u_i_alpha, V_u, @funarg);
k_func = funcreate(@gpc_evaluate, k_i_alpha, V_k, @funarg);

rho=5;
measurements = {
    {@(u)(u(20,:)), gendist_create('normal', {10, rho*0.1}), 'H'};
    {@(u)(u(80,:)), gendist_create('normal', {5, rho*0.3}), 'H'};
    };

[ym_beta, V_ym] = create_measurement_gpc(measurements);
M_func=funcreate(@do_measurement, measurements, @funarg);
y_func = funcompose(u_func, M_func);

[un_i_beta, V_un]=mmse_compose_error_model(u_func, y_func, V_u, ym_beta, V_ym, p_phi, p_int_mmse, p_un, p_int_proj);

n = 20;
mh=multiplot_init(2,1);
multiplot; plot_pce_realizations_1d( pos, u_i_alpha, V_u{2}, 'realizations', n, 'show_stat', 3 );
grid on;
multiplot; plot_pce_realizations_1d( pos, un_i_beta, V_un{2}, 'realizations', n, 'show_stat', 3 );
grid on;
same_scaling(mh, 'y');

save_figure(mh(1), 'u_orig', 'figdir', jb_figdir);
save_figure(mh(2), 'u_updated', 'figdir', jb_figdir);



[kn_i_beta, V_kn]=mmse_compose_error_model(k_func, y_func, V_u, ym_beta, V_ym, p_phi, p_int_mmse, p_un, p_int_proj);

n = 20;
mh=multiplot_init(2,1);
multiplot; plot_pce_realizations_1d( pos, k_i_alpha, V_k{2}, 'realizations', n, 'show_stat', 3 );
grid on;
multiplot; plot_pce_realizations_1d( pos, kn_i_beta, V_kn{2}, 'realizations', n, 'show_stat', 3 );
grid on;
same_scaling(mh, 'y');

save_figure(mh(1), 'k_orig', 'figdir', jb_figdir);
save_figure(mh(2), 'k_updated', 'figdir', jb_figdir);



function [xn_i_beta, V_xn]=mmse_compose_error_model(x_func, y_func, V_xy, ym_beta, V_ym, p_phi, p_int_mmse, p_un, p_int_proj)

[phi_j_delta,V_phi]=mmse_estimate(x_func, y_func, V_xy, p_phi, p_int_mmse);

phi_func = funcreate(@gpc_evaluate, phi_j_delta, V_phi);
ym_func = funcreate(@gpc_evaluate, ym_beta, V_ym);
xn_func = funcompose(ym_func, phi_func);

V_xn = gpcbasis_create(V_ym, 'p', p_un);
xn_i_beta = gpc_projection(xn_func, V_xn, p_int_proj);


function [ym_beta, V_ym] = create_measurement_gpc(Mi)
m = size(Mi,1);
ym_beta = zeros(0,1);
V_ym = gpcbasis_create('');
for i=1:m
    [yi_beta, V_yi] = gpc_param_expand(Mi{i}{2}, Mi{i}{3});
    [ym_beta, V_ym] = gpc_combine_inputs(ym_beta, V_ym, yi_beta, V_yi);
end



function y=do_measurement(Mi, u)
m = size(Mi,1);
n = size(u,2);
y = zeros(m,n);
for i=1:m
    y(i,:) = funcall(Mi{i}{1}, u);
end


function [r_i_alpha, V_r, pos, els] = load_spde_solution(what)
cache_script('sample_solve_spde');
switch what
    case 'u'
        I_r=I_u;
        r_i_alpha = u_i_alpha;
    case 'k'
        I_r=I_k;
        r_i_alpha = kl_to_pce(k_i_k, k_k_alpha);
    otherwise
        error('error:einval', 'unknown: what="%s"', what)
end
V_r = gpcbasis_create('H', 'I', I_r);
pos=pos; %#ok<ASGSL,NODEF>
els=els; %#ok<NODEF,ASGSL>
