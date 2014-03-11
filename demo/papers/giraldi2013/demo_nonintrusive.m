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
N_MC = 1000;
M_MC = 30;
% N_MC = 30;
% M_MC = 3;
params_MC.solve_opts = {'maxiter', 1000, 'steptol', 1e-14, 'abstol', 0};
params_MC.use_surrogate = true;
params_MC.repeat = M_MC;
        
% Define the parameters of the model
p1_dist = gendist_create('uniform', {-1, 1});
p2_dist = gendist_create('uniform', {-1, 1});

% Dishi uses the unnormalised Legendre polynomials (P)
[p1_alpha, V_p1] = gpc_param_expand(p1_dist, 'P');
[p2_alpha, V_p2] = gpc_param_expand(p1_dist, 'P');
[p_alpha, V_p] = gpc_combine_inputs(p1_alpha, V_p1, p2_alpha, V_p2);


results = struct;

for j=1:length(params)
    m = params(j).m;
    underline(strvarexpand('m=$m$'), 'newlines', 1);

    % create basis for polynomial order m
    Vu = gpcbasis_create(V_p, 'p', m);
    %gpcbasis_polynomials(Vu, 'symbols', 'x,y')
    
    % set the parameters according to the paper
    p_int = params(j).p_int;
    steptol = params(j).steptol;
    gal_steptol = params(j).gal_steptol;
    grid = params(j).grid;
    
    % Compute the Collocation response surface and the error
    model = model_stats(model, 'reset');
    [u_alpha_col, model] = compute_response_surface_projection(model, p_alpha, V_p, Vu, p_int, 'steptol', steptol, 'abstol', 0, 'grid', grid);
    model_stats(model, 'print');
    results(j).col_solves = model.model_stats.num_step_calls;
    
    rand('seed', 1234); %#ok<RAND>
    [rmse, rmse_std] = compute_mc_error(fine_model_MC, p_alpha, V_p, u_alpha_col, Vu, N_MC, params_MC);
    results(j).col_error = rmse;
    results(j).col_error_std = rmse_std;

    % Compute the nonint Galerkin response surface and the error
    model = model_stats(model, 'reset');
    [u_alpha_gal, model] = compute_response_surface_nonintrusive_galerkin(model, p_alpha, V_p, Vu, p_int, 'steptol', gal_steptol, 'grid', grid);
    model_stats(model, 'print');
    results(j).gal_solves = model.model_stats.num_step_calls;

    rand('seed', 1234); %#ok<RAND>
    [rmse, rmse_std]= compute_mc_error(fine_model_MC, p_alpha, V_p, u_alpha_gal, Vu, N_MC, params_MC);
    results(j).gal_error = rmse;
    results(j).gal_error_std = rmse_std;

    % Compute the nonint Galerkin response surface and the error
    model = model_stats(model, 'reset');
    [u_alpha_gal, model] = compute_response_surface_nonintrusive_galerkin(model, p_alpha, V_p, Vu, p_int, 'steptol', gal_steptol, 'grid', grid, 'method', 'bfgs');
    model_stats(model, 'print');
    results(j).gab_solves = model.model_stats.num_step_calls;

    rand('seed', 1234); %#ok<RAND>
    [rmse, rmse_std]= compute_mc_error(fine_model_MC, p_alpha, V_p, u_alpha_gal, Vu, N_MC, params_MC);
    results(j).gab_error = rmse;
    results(j).gab_error_std = rmse_std;
end

%% Show the results
underline('Results', 'newlines', 1, 'minwidth', 76);
fprintf('      #iter             error              \n');
fprintf(' m    col  gal  bfgs    colloc            galerkin          gal-bfgs\n');
for j=1:length(results)
    par = params(j);
    res = results(j);
    fprintf( '%2d    %3d  %3d   %3d    %9.3e+-%4.1f%%  %9.3e+-%4.1f%%  %9.3e+-%4.1f%%\n', par.m, ...
        res.col_solves, res.gal_solves, res.gab_solves, ...
        res.col_error, 100*res.col_error_std/res.col_error, ...
        res.gal_error, 100*res.gal_error_std/res.gal_error, ...
        res.gab_error, 100*res.gab_error_std/res.gab_error );
end
