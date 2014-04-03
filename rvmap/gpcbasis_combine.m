function [V, ind_phi1, ind_phi2, ind_xi1, ind_xi2]=gpcbasis_combine(V1, V2, mode)
% GPCBASIS_COMBINE Short description of gpcbasis_combine.
%   GPCBASIS_COMBINE Long description of gpcbasis_combine.
%
% Options
%
% References
%
% Notes
%    'outer' always means that Omega1 and Omega2 are distinct. The germs
%    are distinct and independent. Omega1 x Omega2 is returned.
%    'inner' means that Omega1 and Omega2 are the same. That is the germs
%    are identical.
%    'direct_sum' means 
%
% Example (<a href="matlab:run_example gpcbasis_combine">run</a>)
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

I1 = V1{2};
I2 = V2{2};
switch mode
    case 'outer_sum'
        [I, ind_phi1, ind_phi2, ind_xi1, ind_xi2] = multiindex_outer_direct_sum(I1, I2);
        V = gpcbasis_create(gpcgerm_combine(V1, V2), 'I', I);
    case 'inner_sum'
        [I, ind_phi1, ind_phi2] = multiindex_inner_direct_sum(I1, I2);
        V = gpcbasis_create(V1, 'I', I);
    case 'outer_product'
        [I, ind_phi1, ind_phi2, ind_xi1, ind_xi2] = multiindex_outer_direct_product(I1, I2);
        V = gpcbasis_create(gpcgerm_combine(V1, V2), 'I', I);
    case 'inner_product'
        [I, ind_phi1, ind_phi2] = multiindex_inner_direct_product(I1, I2);
        V = gpcbasis_create(V1, 'I', I);
    otherwise
        error('sglib:gpcbasis_combine', 'unknown combination mode: %s', mode);
end
        

function [I, ind_phi1, ind_phi2, Q1, Q2]=multiindex_outer_direct_sum(I1, I2)
Q1 = 1:size(I1,2);
Q2 = (1:size(I2,2)) + Q1(end);
[I11, I12] = multiindex_combine({I1, I2});
[I, ind_phi1, ind_phi2] = multiindex_inner_direct_sum(I11, I12);

function [I, ind_phi1, ind_phi2]=multiindex_inner_direct_sum(I1, I2)
If = [I1; I2];
I=multiindex_unique(If);
ind_phi1 = multiindex_find(I1, I);
ind_phi2 = multiindex_find(I2, I);

function [I, ind_phi1, ind_phi2, Q1, Q2]=multiindex_outer_direct_product(I1, I2)
Q1 = 1:size(I1,2);
Q2 = (1:size(I2,2)) + Q1;
[I11, I12] = multiindex_combine({I1, I2});
[I, ind_phi1, ind_phi2] = multiindex_inner_direct_product(I11, I12);

function [I, P1, P2]=multiindex_inner_direct_product(I1, I2)
n1 = size(I1,1);
n2 = size(I2,1);
[N1,N2]=meshgrid(1:n1,1:n2);
If = I1(N1(:),:) + I2(N2(:),:);
I=multiindex_unique(If);
P1 = multiindex_find(I1, I);
P2 = multiindex_find(I2, I);

function [I]=multiindex_unique(If)
I=unique(If, 'rows');
I=[sum(I,2), I(:,end:-1:1)];
I=sortrows(I);
I=I(:,end:-1:2);


function g=gpcgerm_combine(V1, V2)
g1 = V1{1};
m1 = gpcbasis_size(V1,2);
g2 = V2{1};
m2 = gpcbasis_size(V2,2);
if all(length(g1)~=[1,m1])
    error('sglib:gpc', 'GPC basis V1 inconsistent');
end
if all(length(g2)~=[1,m2])
    error('sglib:gpc', 'GPC basis V2 inconsistent');
end
if length(g1)==1
    g1 = repmat(g1, 1, m1);
end
if length(g2)==1
    g2 = repmat(g2, 1, m2);
end
g = [g1, g2];
