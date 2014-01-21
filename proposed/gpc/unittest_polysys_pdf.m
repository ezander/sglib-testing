function unittest_polysys_pdf
% UNITTEST_POLYSYS_PDF Test the POLYSYS_PDF function.
%
% Example (<a href="matlab:run_example unittest_polysys_pdf">run</a>)
%   unittest_polysys_pdf
%
% See also POLYSYS_PDF, TESTSUITE 

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

munit_set_function( 'polysys_pdf' );

x=linspace(-3,3);
x=reshape(x, 4, []);
assert_equals( polysys_pdf('h', x), exp(-x.^2/2)/sqrt(2*pi), 'h');
assert_equals( polysys_pdf('p', x), 0.5*(abs(x)<=1), 'p');
assert_equals( polysys_pdf('l', x), (x>=0).*exp(-x), 'l');
assert_equals( polysys_pdf('t', x), 2*(abs(x)<1)./(pi*sqrt(1000*(abs(x)>=1)+(x+1).*(1-x))), 't');
assert_equals( polysys_pdf('u', x), (abs(x)<=1).*(sqrt(1-x.^2))/(pi/4), 'u');

% test the error checking
assert_error( funcreate(@polysys_pdf, 'm',[1,2,3]), 'sglib:gpc', 'no_mono');
assert_error( funcreate(@polysys_pdf, '?',[1,2,3]), 'sglib:gpc', 'unknown');
