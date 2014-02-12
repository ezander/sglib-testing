function sparse_grids_p156(varargin)
%SPARSE_GRIDS_P156 Generates the sparse grid plots from page 156.
%   SPARSE_GRIDS_P156 generates the plots of the sparse grids from [1] page
%   156.
%
% References
%   [1] A. Keese: Numerical Solution of Systems with Stochastic
%       Uncertainties: A General Purpose Framework for Stochastic Finite
%       Elements, Dissertation, TU Braunschweig, 2004,
%       http://digisrv-1.biblio.etc.tu-bs.de:8080/docportal/receive/DocPortal_document_00001595
%
% Example (<a href="matlab:run_example sparse_grids_p156">run</a>)
%    sparse_grids_p156
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

options=varargin2options(varargin,mfilename);
[fejer,options]=get_option(options, 'fejer', false);
check_unsupported_options(options);

multiplot_init(3,4);

rule_func = @gauss_legendre_rule;
rule_name = 'Gauss-Legendre';
show_quad_grid(2, 3, rule_func, rule_name)
show_quad_grid(2, 5, rule_func, rule_name)
show_quad_grid(2, 6, rule_func, rule_name)

rule_func = @gauss_hermite_rule;
rule_name = 'Gauss-Hermite';
show_quad_grid(2, 3, rule_func, rule_name)
show_quad_grid(2, 5, rule_func, rule_name)
show_quad_grid(2, 6, rule_func, rule_name)

if ~fejer
    rule_func = @clenshaw_curtis_nested;
    rule_name = 'Clenshaw-Curtis';
    k = [4,7,9];
else
    rule_func = @(n)(clenshaw_curtis_nested( n, 'mode', 'fejer2'));
    rule_name = 'Fejer 2';
    k = [4,5,7];
end
show_quad_grid(2, k(1), rule_func, rule_name)
show_quad_grid(2, k(2), rule_func, rule_name)
show_quad_grid(2, k(3), rule_func, rule_name)

show_quad_grid(3, k(1), rule_func, rule_name)
show_quad_grid(3, k(2), rule_func, rule_name)
show_quad_grid(3, k(3), rule_func, rule_name)

function show_quad_grid(d, stages, rule_func, rule_name) %#ok<INUSD>
multiplot; 
[x,w] = smolyak_grid(d, stages, rule_func); %#ok<NASGU>
strvarexpand('$rule_name$, d=$d$, k=$stages$: n=$length(w)$')
if d==2
    plot(x(1,:), x(2,:), '.');
    axis square;
else
    plot3(x(1,:), x(2,:), x(3,:), '.');
    axis square;
    view(3);
end
titlestr = strvarexpand('$rule_name$, S^$d$_$stages$, $size(x,2)$ nodes.');
title(titlestr);
