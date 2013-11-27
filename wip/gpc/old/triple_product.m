function [M,p]=triple_product( p, raw_moments_func )
% TRIPLE_PRODUCT Compute triple products for polynomials system.
%     TRIPLE_PRODUCT( P, RAW_MOMENTS_FUNC ) Computes the expectation of the
%     triple products E[p_i,p_j,p_k] for the polynomials given in P. If you
%     pass an integer in P the orthogonal polynomial system up to degree P
%     is computed by the RAW_MOMENTS_FUNC and returned in the output
%     parameter P.
% 
% Example (<a href="matlab:run_example triple_product">run</a>)
%     moments_func = {@normal_raw_moments, {0, 1}, {2, 3}};
%     save_format( 'compact', 'short g' );
%     M=triple_product( 5, moments_func );
%     chopabs(M)
%     restore_format();
%
% See also 

%   Elmar Zander
%   Copyright 2009, Institute of Scientific Computing, TU Braunschweig.
%   $Id$ 
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

if isscalar(p)
    n=p;
    p=orthpoly_from_moments( raw_moments_func, n );
else
    n=size(p,1)-1;
end
m=funcall( raw_moments_func,  0:3*n );
    
M=zeros(n+1,n+1,n+1);
for i=1:n+1
    for j=1:n+1
        for k=1:n+1
            p_ijk=conv( conv( p(i,:), p(j,:)), p(k,:));
            M(i,j,k)=sum( m.*fliplr( p_ijk ) );
        end
    end
end
