function unittest_gpcgerm_pdf
% UNITTEST_GPCGERM_PDF Test the GPCGERM_PDF function.
%
% Example (<a href="matlab:run_example unittest_gpcgerm_pdf">run</a>)
%   unittest_gpcgerm_pdf
%
% See also GPCGERM_PDF, TESTSUITE 

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

munit_set_function( 'gpcgerm_pdf' );

x = linspace(-3,3);
y = rand(size(x));

V=gpcbasis_create('h');
assert_equals(gpcgerm_pdf(V,x), polysys_pdf('h', x), 'h');

V=gpcbasis_create('ul');
assert_equals(gpcgerm_pdf(V,[x;y]), polysys_pdf('u', x).*polysys_pdf('l', y), 'ul');

V=gpcbasis_create('h', 'm', 2);
assert_equals(gpcgerm_pdf(V,[x;y]), polysys_pdf('h', x).*polysys_pdf('h', y), 'h2');

