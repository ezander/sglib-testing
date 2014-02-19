function block_sparsity_p72(varargin)
%BLOCK_SPARSITY_P72 Generates the sparsity plots from page 72.
%   BLOCK_SPARSITY_P72 Generates the sparsity plots from A. Keese's diss.
%   (fig 41. on page 72). 
%
% Example (<a href="matlab:run_example block_sparsity_p72">run</a>)
%   block_sparsity_p72
%
% See also MULTIINDEX, HERMITE_TRIPLE_FAST

%   Elmar Zander
%   Copyright 2009, Institute of Scientific Computing, TU Braunschweig.
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

options=varargin2options(varargin,mfilename);
[show_other_ordering,options]=get_option(options, 'show_other_ordering', false);
check_unsupported_options(options);

sp_plots('lex', false, 'Sparsity plot (lexicographic ordering, A. Keese)');

if show_other_ordering
    userwait;
    sp_plots('degree', false, 'Sparsity plot (degree ordering)');

    userwait;
    sp_plots('uqtk', false, 'Sparsity plot (uq toolkit)');

    userwait;
    sp_plots('degree', true, 'Sparsity plot (degree + reduced bandwidth ordering)');
end



function sp_plots(ordering, do_rcm, titlestr )
% SP_PLOTS Show sparsity plots.
% If DO_SORTING is true then the multiindices are sorted lexicographically
% like in A. Keese's thesis. Otherwise a (to me more natural) ordering is
% used, where lower orders come first, and for each order, the multiindices
% are then ordered lexicographically (highest ranking right). If DO_RCM is
% true a bandwidth reduction based on the reverse Cuthill-McKee algorithm
% is used (DO_SORTING is more or less no effect then).

clf;
m=4;
p=4;
p2=[1,4,2,5];
I=multiindex(m,p,'ordering', ordering);

V1 = gpcbasis_create('H', 'I', I);
for i=1:4
    J=multiindex(m, p2(i), 'ordering', ordering);
    V2 = gpcbasis_create('H', 'I', J);
    
    % compute sparsity pattern
    C = gpcbasis_triples(V1, V1, V2);
    S = any(C, 3);
    subplot(2,2,i);

    % bandwidth reduction with reverse cuthill-mckee
    if do_rcm
        r=symrcm(S);
        S=S(r,r);
    end;
    
    spy(S);
    drawnow;
    
    %print bandwidth
    [i1,i2]=find(sparse(S));
    bw = sprintf( 'bandwidth: %d\n', 1+max(abs(i1-i2)) ); %#ok<NASGU>
    title(strvarexpand('p: $p2(i)$, $bw$'));
end
annotation('textbox', [0 0.89 1 0.1], ...
    'String', titlestr, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
