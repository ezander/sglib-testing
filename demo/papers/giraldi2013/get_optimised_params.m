function params = get_optimised_params()
% Set the parameters according to the paper 

params = struct();
for j = 1:4
    % degree goes from 2 to 5
    m = j+1;
    params(j).m = m;
    
    % integration order one more than polynomial order    
    params(j).p_int = m+1;
    
    % step tolerance for the solvers (this seems to be enough)
    params(j).steptol = 10^-(m+3);
    
    % for the stochastic Galerkin solver this needs to be modified (this is
    % not from the paper, but from the actual code, because here the
    % function mean() is used instead of norm())
    n_space = 5;
    n_stoch = multiindex_size(2, m);
    n_dofs = n_space * n_stoch;
    params(j).gal_steptol = params(j).steptol * sqrt(n_dofs);
    
    % grid is always a full tensor grid
    params(j).grid = 'full_tensor';
end
