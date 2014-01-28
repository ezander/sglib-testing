function A_i = gpc_multiplication_matrices(a_i_alpha, V_a, V_u)
% GPC_MULTIPLICATION_MATRICES Short description of gpc_multiplication_matrices.
%   GPC_MULTIPLICATION_MATRICES Long description of gpc_multiplication_matrices.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example gpc_multiplication_matrices">run</a>)
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

M_a = gpcbasis_size(V_a, 1);
M_u = gpcbasis_size(V_u, 1);

M = reshape(gpc_triples(V_a, V_u, V_u), M_a, []);
A = a_i_alpha * M;
R = size(a_i_alpha,1);
A_i = cell(1,R);
for i=1:R
    A_i{i} = reshape(A(i,:), M_u, M_u);
end

