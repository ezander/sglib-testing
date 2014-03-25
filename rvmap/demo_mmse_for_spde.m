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

[u_i_alpha, V_u, pos, els] = load_pde_solution;

x_func = funcreate(@gpc_evaluate, u_i_alpha, V_u, @funarg);
y_func = funcompose(x_func, @measurement);

[phi_j_delta,V_phi]=mmse_estimate_generic(x_func, y_func, V_u, 3, 3);


ym1_dist = gendist_create('normal', {5,2});
ym2_dist = gendist_create('normal', {15,1});
ym1_dist = gendist_create('normal', {5,3.2});
ym2_dist = gendist_create('normal', {13,0.3});


ym1_dist = gendist_create('normal', {1.3,0.01});
ym2_dist = gendist_create('normal', {0.2,0.01});


[ym1_beta, V_m1] = gpc_param_expand(ym1_dist, 'H');
[ym2_beta, V_m2] = gpc_param_expand(ym2_dist, 'H');
[ym_beta, V_ym] = gpc_combine_inputs(ym1_beta, V_m1, ym2_beta, V_m2);


phi_func = funcreate(@gpc_evaluate, phi_j_delta, V_phi);
ym_func = funcreate(@gpc_evaluate, ym_beta, V_ym);
phi_ym_func = funcompose(ym_func, phi_func);

p_un = 3;
p_int = 6;
[un_i_beta, V_un] = compute_gpc_projection(phi_ym_func, V_ym, p_int, p_un)

n = 20;
multiplot; plot_pce_realizations_1d( pos, un_i_beta, V_un{2}, 'realizations', n, 'show_stat', 3 );
grid on;


function func=funcompose(func1, func2)
func = funcreate(@composed_function, func1, func2, @funarg);

function z=composed_function(func1, func2, x)
y = funcall(func1, x);
z = funcall(func2, y);


function y=measurement(x)
%y = [x(6, :); x(77,:)];
y = [x(12, :); x(77,:)];

function [un_i_beta, V_un] = compute_gpc_projection(un_func, V_u, p_int, p_un)
[xi_k,w_k] = gpc_integrate([], V_u, p_int);
Q = length(w_k);

V_un = gpcbasis_create(V_u, 'p', p_un);
M = gpcbasis_size(V_un, 1);
un_i_beta = zeros(100, M);

for k = 1:Q
    u_i_j = funcall(un_func, xi_k(:,k));
    psi_j_beta_dual = gpcbasis_evaluate(V_un, xi_k(:,k), 'dual', true);
    un_i_beta = un_i_beta + w_k(k) * u_i_j * psi_j_beta_dual;
end





function [r_i_alpha, V_r, pos, els] = load_pde_solution
cache_script('sample_solve_spde');
I_r=I_u;
r_i_alpha = u_i_alpha;
I_r=I_k;
r_i_alpha = kl_to_pce(k_i_k, k_k_alpha);
V_r = gpcbasis_create('H', 'I', I_r);
pos;
els;

% plot solutions fields and difference
mh=multiplot_init(2,1);
switch(d)
    case 1
        n = 20;
        multiplot; 
        plot_pce_realizations_1d( pos, r_i_alpha, I_r, 'realizations', n, 'show_stat', 3 ); grid on;
    case 2
        opts={'view', 3};
        multiplot; plot_field(pos, els, r, opts{:} );
end




function [u_i_alpha, V_u, pos, els] = load_pde_solution_old
cache_script('sample_solve_spde');
V_u = gpcbasis_create('H', 'I', I_u);
u_i_alpha;
pos;
els;

% plot solutions fields and difference
mh=multiplot_init(2,2);
switch(d)
    case 1
        n = 20;
        multiplot; plot_kl_pce_realizations_1d( pos, k_i_k, k_k_alpha, I_k, 'realizations', n, 'show_stat', 3 );
        multiplot; plot_kl_pce_realizations_1d( pos, f_i_k, f_k_alpha, I_f, 'realizations', n, 'show_stat', 3 );
        multiplot; plot_pce_realizations_1d( pos, u_i_alpha, I_u, 'realizations', n, 'show_stat', 3 ); grid on;
    case 2
        opts={'view', 3};
        multiplot; plot_field(pos, els, u, opts{:} );
        multiplot; plot_field(pos, els, k, opts{:} );
        multiplot; plot_field(pos, els, f, opts{:} );
        multiplot; plot_field(pos, els, g, opts{:} );
end

