function unittest_orthpoly_from_moments
% UNITTEST_ORTHPOLY_FROM_MOMENTS Test the ORTHPOLY_FROM_MOMENTS function.
%
% Example (<a href="matlab:run_example unittest_orthpoly_from_moments">run</a>)
%   unittest_orthpoly_from_moments
%
% See also ORTHPOLY_FROM_MOMENTS 

%   Elmar Zander
%   Copyright 2009, Inst. of Scientific Computing, TU Braunschweig
%   $Id$ 
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.


munit_set_function( 'orthpoly_from_moments' );


expect=[
    0    0    0    0    0    1
    0    0    0    0    1    0
    0    0    0    1    0   -1
    0    0    1    0   -3    0
    0    1    0   -6    0    3
    1    0  -10    0   15    0];
p=orthpoly_from_moments( {@normal_raw_moments, {0, 1}, {2, 3}}, 5, 'monic', true );
assert_equals( p, expect, 'H_monic' );


expect=[
    0, 0, 0, 0, 1
    0, 0, 0, 1, 0
    0, 0, 1, 0, -1/3
    0, 1, 0, -3/5, 0
    1, 0, -6/7, 0, 3/35];

p=orthpoly_from_moments( {@uniform_raw_moments, {-1, 1}, {2, 3}}, 4, 'monic', true );
assert_equals( p, expect, 'L_monic' );



%p=orthpoly_from_moments( {@lognormal_raw_moments, {0, 1}, {1, 2}}, 5 );
