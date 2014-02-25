clear

% In contrast to what the paper says R is 1/100 not 100
model = electrical_network_init('R', 1/100, 'newton', false);

a1_dist = gendist_create('uniform', {-1, 1});
a2_dist = gendist_create('uniform', {-1, 1});

% Dishi uses the unnormalised Legendre polynomials
[a1_alpha, Va1] = gpc_param_expand(a1_dist, 'P');
[a2_alpha, Va2] = gpc_param_expand(a1_dist, 'P');
[a_alpha, Va] = gpc_combine_inputs(a1_alpha, Va1, a2_alpha, Va2);


clc
N_mc = 100;
fine_model = electrical_network_init('newton', true);
results = struct;
for j=1:4
    % create basis for polynomial order m = j+1
    m = j + 1;
    results(j).m = m;
    Vu = gpcbasis_create(Va, 'p', m);
    %gpcbasis_polynomials(Vu, 'symbols', 'x,y')
    
    % set the parameters according to the paper: integration order one more
    % than polynomial order, the tolerance needs to be scaled by sqrt(n),
    % i.e. the total degrees of freedom (this is not from the paper, but
    % from the actual code), and set the grid to full tensor
    p_int = m+1;
    steptol = 10^-(m+4);
    n = model.model_info.num_vars * gpcbasis_info(Vu, 'num_basis_functions');
    gal_steptol = steptol * sqrt(n);
    grid = 'full_tensor';
    
    % Now compute the nonint Galerkin response surface and the error
    underline(strvarexpand('m=$m$'), 'newlines', 1);
    model = model_stats(model, 'reset');
    [u_alpha_gal, model] = compute_response_surface_nonintrusive_galerkin(model, a_alpha, Va, Vu, p_int, 'steptol', gal_steptol, 'grid', grid);
    model_stats(model, 'print');
    results(j).gal_solves = model.model_stats.num_step_calls;

    rand('seed', 1234); %#ok<RAND>
    rmse = compute_mc_error(fine_model, a_alpha, Va, u_alpha_gal, Vu, N_mc, 'solve_opts', {'maxiter', 1000, 'steptol', 1e-14, 'abstol', 0});
    results(j).gal_error = rmse;
    
    % Then compute the Collocation response surface and the error
    model = model_stats(model, 'reset');
    [u_alpha_col, model] = compute_response_surface_projection(model, a_alpha, Va, Vu, p_int, 'steptol', steptol, 'abstol', 0, 'grid', grid);
    model_stats(model, 'print');
    results(j).col_solves = model.model_stats.num_step_calls;
    
    rand('seed', 1234); %#ok<RAND>
    rmse = compute_mc_error(fine_model, a_alpha, Va, u_alpha_col, Vu, N_mc, 'solve_opts', {'maxiter', 1000, 'steptol', 1e-14, 'abstol', 0});
    results(j).col_error = rmse;
end

% Show the results
underline('Results', 'newlines', 1);
strvarexpand('m    col    gal    col           gal       ');
for j=1:length(results)
    r = results(j);
    strvarexpand('$r.m$    $r.col_solves$    $r.gal_solves$    $r.col_error$    $r.gal_error$');
end
