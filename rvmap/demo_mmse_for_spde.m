function demo_mmse_for_spde(varargin)
% DEMO_MMSE_FOR_SPDE Demonstration of the MMSE methods for a sample SPDE.
%
% See also MMSE_ESTIMATE, MMSE_UPDATE_MODEL

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

% Set up the parameters for integration order and expansion degrees
p_phi = 3;
p_int_mmse = 4;
p_un = 4;
p_int_proj = p_un+1;

% Load the SPDE surrogate models for U and KAPPA and construct evaluation
% functions for them
[u_i_alpha, V_u, pos] = load_spde_solution('u');
[k_i_alpha, V_k] = load_spde_solution('k');

u_func = gpc_function(u_i_alpha, V_u);
k_func = gpc_function(k_i_alpha, V_k);

% Construct a struct for the first measurement (measuring u at 0.2 giving a
% value of u(0.2)=10+-0.5 (normally distributed)
M(1).func = @(u)(u(30,:));
M(1).dist = gendist_create('normal', {22, 0.2});
M(1).gpc_polysys = 'H';

% Construct a struct for the second measurement (measuring u at 0.8 giving a
% value of u(0.2)=5+-0.5 (normally distributed);
M(2).func = @(u)(u(50,:));
M(2).dist = gendist_create('normal', {28, 10.3});
M(2).gpc_polysys = 'H';

% Construct a struct for the second measurement (measuring u at 0.8 giving a
% value of u(0.2)=5+-0.5 (normally distributed);
M(3).func = @(u)(u(80,:));
M(3).dist = gendist_create('normal', {18, 0.3});
M(3).gpc_polysys = 'H';

% Create a function for the measurement operator M
M_func=funcreate(@do_measurement, {M.func}, @funarg);

% Build the composition y = M o u
y_func = funcompose(u_func, M_func);

% Update the stochastic model for U using the measurements y_m
% Constructs a GPC of the combined actual measurements plus error model y_m
[ym_beta, V_ym] = create_measurement_gpc({M.dist}, {M.gpc_polysys});
ym = ym_beta(:,1);
eps_j_gamma = ym_beta;
eps_j_gamma(:,1) = 0;
V_eps = V_ym;
eps_func = gpc_function(eps_j_gamma, V_eps);
p_pn = p_un;


[un_i_beta, V_un]=mmse_update_gpc(  u_i_alpha, y_func, V_u, ym, eps_func, V_eps, p_phi, p_int_mmse, p_pn, p_int_proj);




% Plot the original and updated U
n = 20;
mh=multiplot_init(2,1);
multiplot; plot_pce_realizations_1d( pos, u_i_alpha, V_u{2}, 'realizations', n, 'show_stat', 3 );
grid on;
multiplot; plot_pce_realizations_1d( pos, un_i_beta, V_un{2}, 'realizations', n, 'show_stat', 3 );
grid on;
same_scaling(mh, 'y');

save_figure(mh(1), 'u_orig', 'figdir', jb_figdir);
save_figure(mh(2), 'u_updated', 'figdir', jb_figdir);


% Update the stochastic model for KAPPA using the measurements y_m
[kn_i_beta, V_kn]=mmse_update_gpc(  k_i_alpha, y_func, V_u, ym, eps_func, V_eps, p_phi, p_int_mmse, p_pn, p_int_proj);

% Plot the original and updated KAPPA
n = 20;
mh=multiplot_init(2,1);
multiplot; plot_pce_realizations_1d( pos, k_i_alpha, V_k{2}, 'realizations', n, 'show_stat', 3 );
grid on;
multiplot; plot_pce_realizations_1d( pos, kn_i_beta, V_kn{2}, 'realizations', n, 'show_stat', 3 );
grid on;
same_scaling(mh, 'y');

save_figure(mh(1), 'k_orig', 'figdir', jb_figdir);
save_figure(mh(2), 'k_updated', 'figdir', jb_figdir);



function [r_i_alpha, V_r, pos, els] = load_spde_solution(what)
[u_i_alpha, I_u, k_i_alpha, I_k, pos, els] = load_complete_spde_solution();
switch what
    case 'u'
        I_r=I_u;
        r_i_alpha = u_i_alpha;
    case 'k'
        I_r=I_k;
        r_i_alpha = k_i_alpha;
    otherwise
        error('error:einval', 'unknown: what="%s"', what)
end
V_r = gpcbasis_create('H', 'I', I_r);

function [u_i_alpha, I_u, k_i_alpha, I_k, pos, els] = load_complete_spde_solution() %#ok<STOUT>
cache_script('sample_solve_spde');
k_i_alpha = kl_to_pce(k_i_k, k_k_alpha);
