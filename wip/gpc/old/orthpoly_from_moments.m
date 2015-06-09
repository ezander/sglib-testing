function p=orthpoly_from_moments( raw_moments_func, n, varargin )
% ORTHPOLY_FROM_MOMENTS Compute system of orthogonal polynomials.
%   ORTHPOLY_FROM_MOMENTS( RAW_MOMENTS_FUNC, N ) computes the system of
%   orthogonal polynomials up to degree N (i.e. N+1 polynomials in sum) for
%   the distribution specified throught its raw moments (given in
%   RAW_MOMENTS_FUNC).
%
% Options:
%   monic: true, {false}
%     Returns the polynomials not normalized (which is the default) but as
%     monic polynomials with leading coefficient 1. This makes it easier to
%     compare to tables of polynomials in reference works.
% 
% Example 1 (<a href="matlab:run_example orthpoly_from_moments">run</a>)
%     underline( 'The first 6 Hermite polynomials (up to degree 5)' );
%     p=orthpoly_from_moments( {@normal_raw_moments, {0, 1}, {2, 3}}, ...
%       5, 'monic', true );
%     format_poly(p,'tight',false,'twoline',true,'symbol','t')
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

options=varargin2options( varargin );
[monic,options]=get_option( options, 'monic', false );
check_unsupported_options( options );


% Description of the algorithm:
% Some easy calculation shows that if 
% $\int p_i p_j d\mu(x)=delta_{i,j}$ then the coefficient vectors $a_i$ and
% $a_j$ are orthogonal with respect to the matrix 
% $M_{i,j}=\int x^{i+j} d\mu(x)=m_{i+j}$ (zero based) where $m_k$ denotes
% the k-th moment of the distribution given by $\mu$. Thus we can create
% $M$ as a Hankel matrix of the raw moments of the distribution, then
% orthogonalize the monomials (represented as the identity matrix) with
% respect to that matrix.

% create the moment matrix, 
raw_moments=funcall( raw_moments_func, 0:2*n );
M=hankel( raw_moments(1:n+1), raw_moments(n+1:end) );

% the following gets us the orthogonal polynomials as column vectors
% p=gram_schmidt( eye(n+1), M );
p = inv(chol(M, 'upper'));

% however, we need them as row vectors and in the braindead matlab right to
% left ordering
p=fliplr( p');

% if monic polynomials are requested, divide by the leading term
if monic
    p=diag(diag(fliplr(p)))\p;
end

