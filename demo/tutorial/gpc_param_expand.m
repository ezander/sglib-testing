function [a_alpha, V] = gpc_param_expand(dist_a, sys, p,varargin)
% GPC_PARAM_EXPAND Short description of gpc_param_expand.
%   GPC_PARAM_EXPAND Long description of gpc_param_expand.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example gpc_param_expand">run</a>)
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

p_int = max(12, 2*p);

V = gpcbasis_create(sys, 'p', p);

[x,w]=gpc_integrate([], V, p_int);

psi_k_alpha = gpcbasis_evaluate(V, x, 'dual', true);
fun_k = gendist_invcdf(gpcgerm_cdf(V, x), dist_a{:});
a_alpha = fun_k*diag(w)*psi_k_alpha;
