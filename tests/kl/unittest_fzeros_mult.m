function unittest_fzeros_mult
% UNITTEST_FZEROS_MULT Test the FZEROS_MULT function.
%
% Example (<a href="matlab:run_example unittest_fzeros_mult">run</a>)
%   unittest_fzeros_mult
%
% See also FZEROS_MULT, MUNIT_RUN_TESTSUITE 

%   Elmar Zander
%   Copyright 2013, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

munit_set_function( 'fzeros_mult' );

% detect all zeros of a polynomial
r0 = [0.1,0.2,0.1*pi];
func=funcreate(@polyval, poly(r0), @funarg);
r=fzeros_mult( func, 0, 10 );
assert_equals(r, r0, 'poly');

% its ok to not detect double zeros
r0 = [0.1,0.1,0.2,0.3];
func=funcreate(@polyval, poly(r0), @funarg);
r=fzeros_mult( func, 0, 10, 'N', 1000 );
assert_equals(r, r0(3:end), 'poly_sign_change');

% don't report minima and sign changes at poles
z=0.1*pi;
func=@(x)( ((x-z).^2+1) .* (x-2*z) .* (x-4*z) ./ (x-3*z) );
r=fzeros_mult( func, 0, 10 );
assert_equals(r, [2*z, 4*z], 'min_and_inf');

% a special function (there should be N zeros in [0,N*pi])
fun1=@(x)(1-x.*tan(x));
r=fzeros_mult( fun1, 0, 7*pi, 'N', 100 );
assert_equals(size(r), [1,7], 'num');
assert_equals(fun1(r), zeros(size(r)), 'zero_tan');

% same special function (compare to approximation)
r=fzeros_mult( fun1, 90*pi-0.1, 100*pi );
N=90:99;
assert_equals(r, pi*(0.0011 + N), 'tan_approx', 'reltol', 1e-6);
