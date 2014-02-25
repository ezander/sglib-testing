% DEMO_NONINSTRUSIVE Recreate tables published in non-intrusive Galerkin paper.
% This demo tries to reproduce the results from [1]
%
% [1] L. Giraldi, A. Litvinenko, D. Liu, H.G. Matthies, A. Nouy: To be or
%     not to be intrusive? The solution of parametric and stochastic
%     equations - the “plain vanilla” Galerkin case Informatikbericht Nr.
%     2013-03, Inst. of Scientific Computing, TU Braunschweig

clear; clc;

% Get simulation parameters
params = get_paper_params();

% In contrast to what the paper says R is 1/100 not 100. If yo set
% use_newton to true, iterations will a) be faster b) converge also for
% higher values of R.
R = 1/100;
fg = 25;
use_newton = false;
model = electrical_network_init('R', R, 'fg', fg, 'newton', use_newton);

% Model and number of evaluations for MC for error estimation (MC may
% always use Newton) (I take 100 here instead of 1000, as that makes it 10
% times faster and there it doesn't seem to make a big difference).
fine_model_MC = electrical_network_init('R', R, 'fg', fg, 'newton', true);
N_MC = 100;
params_MC.solve_opts = {'maxiter', 1000, 'steptol', 1e-14, 'abstol', 0};
params_MC.use_surrogate = true;
        
% Define the parameters of the model
p1_dist = gendist_create('uniform', {-1, 1});
p2_dist = gendist_create('uniform', {-1, 1});

% Dishi uses the unnormalised Legendre polynomials (P)
[p1_alpha, V_p1] = gpc_param_expand(p1_dist, 'P');
[p2_alpha, V_p2] = gpc_param_expand(p1_dist, 'P');
[p_alpha, V_p] = gpc_combine_inputs(p1_alpha, V_p1, p2_alpha, V_p2);


results = struct;

for j=1:length(params)
    % create basis for polynomial order m
    m = params(j).m;
    Vu = gpcbasis_create(V_p, 'p', m);
    %gpcbasis_polynomials(Vu, 'symbols', 'x,y')
    
    % set the parameters according to the paper
    p_int = params(j).p_int;
    steptol = params(j).steptol;
    gal_steptol = params(j).gal_steptol;
    grid = params(j).grid;
    
    % Now compute the nonint Galerkin response surface and the error
    underline(strvarexpand('m=$m$'), 'newlines', 1);
    model = model_stats(model, 'reset');
    [u_alpha_gal, model] = compute_response_surface_nonintrusive_galerkin(model, p_alpha, V_p, Vu, p_int, 'steptol', gal_steptol, 'grid', grid);
    model_stats(model, 'print');
    results(j).gal_solves = model.model_stats.num_step_calls;

    rand('seed', 1234); %#ok<RAND>
    rmse = compute_mc_error(fine_model_MC, p_alpha, V_p, u_alpha_gal, Vu, N_MC, params_MC);
    results(j).gal_error = rmse;
    
    % Then compute the Collocation response surface and the error
    model = model_stats(model, 'reset');
    [u_alpha_col, model] = compute_response_surface_projection(model, p_alpha, V_p, Vu, p_int, 'steptol', steptol, 'abstol', 0, 'grid', grid);
    model_stats(model, 'print');
    results(j).col_solves = model.model_stats.num_step_calls;
    
    rand('seed', 1234); %#ok<RAND>
    rmse = compute_mc_error(fine_model_MC, p_alpha, V_p, u_alpha_col, Vu, N_MC, params_MC);
    results(j).col_error = rmse;
end

% Show the results
underline('Results', 'newlines', 1);
fprintf('    #iter            error              \n');
fprintf(' m  col  gal   colloc     galerkin       \n');
for j=1:length(results)
    par = params(j);
    res = results(j);
    fprintf( '%2d  %3d  %3d  %9.3e  %9.3e\n', par.m, res.col_solves, ...
        res.gal_solves, res.col_error, res.gal_error );
end
