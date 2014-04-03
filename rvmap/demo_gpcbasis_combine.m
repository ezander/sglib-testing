function demo_gpcbasis_combine(varargin)
% DEMO_GPCBASIS_COMBINE Short description of demo_gpcbasis_combine.
%   DEMO_GPCBASIS_COMBINE Long description of demo_gpcbasis_combine.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example demo_gpcbasis_combine">run</a>)
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

clc

% Outer sum and product
V1 = gpcbasis_create('MM', 'p', 2);
V2 = gpcbasis_create('M', 'p', 3);

basis1 = gpcbasis_polynomials(V1, 'symbols', 'xy')'
basis2 = gpcbasis_polynomials(V2, 'symbols', 'z')'

V = gpcbasis_combine(V1, V2, 'outer_sum');
outer_sum = gpcbasis_polynomials(V, 'symbols', 'xyz')'

V = gpcbasis_combine(V1, V2, 'outer_product');
outer_product = gpcbasis_polynomials(V, 'symbols', 'xyz')'

% Inner sum and product
V1 = gpcbasis_create('MM', 'p', 2);
V2 = gpcbasis_create('MM', 'I', [0, 0; 1, 4]);

basis1 = gpcbasis_polynomials(V1, 'symbols', 'xy')'
basis2 = gpcbasis_polynomials(V2, 'symbols', 'xy')'

V = gpcbasis_combine(V1, V2, 'inner_sum');
inner_sum = gpcbasis_polynomials(V, 'symbols', 'xy')'

V = gpcbasis_combine(V1, V2, 'inner_product');
inner_product = gpcbasis_polynomials(V, 'symbols', 'xy')'

