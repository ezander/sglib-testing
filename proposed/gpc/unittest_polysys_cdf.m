function unittest_polysys_cdf
% UNITTEST_POLYSYS_CDF Test the POLYSYS_CDF function.
%
% Example (<a href="matlab:run_example unittest_polysys_cdf">run</a>)
%   unittest_polysys_cdf
%
% See also POLYSYS_CDF, TESTSUITE 

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

munit_set_function( 'polysys_cdf' );

x=linspace(-3,3);
x=reshape(x, 4, []);
assert_equals( polysys_cdf('h', x), 0.5*(1+erf(x/sqrt(2))), 'h');
assert_equals( polysys_cdf('p', x), min(1,max(0,0.5*(1+x))), 'p');
assert_equals( polysys_cdf('l', x), (x>=0).*(1-exp(-x)), 'l');


x=linspace(-0.9,0.9,8);
x=reshape(x, 4, []);
assert_equals( polysys_cdf('t', x), (abs(x)<1).*(asin(x)/pi+0.5), 't');
%assert_equals( polysys_cdf('u', x), (abs(x)<=1).*(sqrt(1-x.^2))/(pi/4), 'u');

% test the error checking
assert_error( funcreate(@polysys_cdf, 'm',[1,2,3]), 'sglib:gpc', 'no_mono');
assert_error( funcreate(@polysys_cdf, '?',[1,2,3]), 'sglib:gpc', 'unknown');
