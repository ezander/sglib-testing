function u_i_beta = gpc_projection(u_func, V_u, p_int)
% GPC_PROJECTION Short description of gpc_projection.
%   GPC_PROJECTION Long description of gpc_projection.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example gpc_projection">run</a>)
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

[xi_k,w_k] = gpc_integrate([], V_u, p_int);
Q = length(w_k);

M = gpcbasis_size(V_u, 1);
u_i_beta = zeros(100, M);

for k = 1:Q
    u_i_j = funcall(u_func, xi_k(:,k));
    psi_j_beta_dual = gpcbasis_evaluate(V_u, xi_k(:,k), 'dual', true);
    u_i_beta = u_i_beta + w_k(k) * u_i_j * psi_j_beta_dual;
end
