function unittest_myfft(varargin)
% UNITTEST_MYFFT Test the MYFFT function.
%
% Example (<a href="matlab:run_example unittest_myfft">run</a>)
%   unittest_myfft
%
% See also MYFFT, MUNIT_RUN_TESTSUITE 

%   Elmar Zander
%   Copyright 2015, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

munit_set_function( 'myfft' );


func = @(x)(exp(x));
[A_k, wp_k, x_i] = myfft(func, -1, 1, 10);
assert_equals(sin_eval(A_k, wp_k, x_i), func(x_i))

[A_k, wp_k, x_i] = myfft(func, -1, 1, 11);
assert_equals(sin_eval(A_k, wp_k, x_i), func(x_i))

[A_k, wp_k, x_i] = myfft(func, -2, 2, 11);
assert_equals(sin_eval(A_k, wp_k, x_i), func(x_i))

[A_k, wp_k, x_i] = myfft(func, 0, 2, 11);
assert_equals(sin_eval(A_k, wp_k, x_i), func(x_i))

[A_k, wp_k, x_i] = myfft(func, 0.01, 2.5, 141);
assert_equals(sin_eval(A_k, wp_k, x_i), func(x_i))
%plot(x_i, func(x_i), x_i, sin_eval(A_k, wp_k, x_i))



func = @(x)(exp(-abs(x)));
[A_k, wp_k, x_i] = myfft(func, -1, 1, 20, 'symmetry', 'even');
assert_equals(sin_eval(A_k, wp_k, x_i), func(x_i))


func = @(x)(sign(x).*(exp(-abs(x))-1));
[A_k, wp_k, x_i] = myfft(func, -1, 1, 100, 'symmetry', 'odd');
x_i = x_i(2:end);
assert_equals(sin_eval(A_k, wp_k, x_i), func(x_i))
