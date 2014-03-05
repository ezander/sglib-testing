function unittest_BetaDistribution
% UNITTEST_BETADISTRIBUTION Test the BETADISTRIBUTION function.
%
% Example (<a href="matlab:run_example unittest_BetaDistribution">run</a>)
%   unittest_BetaDistribution
%
% See also BETADISTRIBUTION, MUNIT_RUN_TESTSUITE 

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

munit_set_function( 'BetaDistribution' );


dist = BetaDistribution(3, 5);
assert_equals(dist.a, 3, 'init a');
assert_equals(dist.b, 5, 'init b');


dist = BetaDistribution(2, 6);
x = rand(5);
assert_equals(dist.pdf(x), beta_pdf(x,2,6), 'pdf');

dist = BetaDistribution(2, 2);
assert_equals(dist.pdf(0.5), 3/2, 'pdf');


