function unittest_fourier_series_expand(varargin)
% UNITTEST_FOURIER_SERIES_EXPAND Test the FOURIER_SERIES_EXPAND function.
%
% Example (<a href="matlab:run_example unittest_fourier_series_expand">run</a>)
%   unittest_fourier_series_expand
%
% See also FOURIER_SERIES_EXPAND, MUNIT_RUN_TESTSUITE 

%   Elmar Zander
%   Copyright 2015, Inst. of Scientific Computing
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

munit_set_function( 'fourier_series_expand' );

func = @(x)(exp(x));
[A_k, TB, x_i] = fourier_series_expand(func, -1, 1, 10);
assert_equals(trig_eval(A_k, TB, x_i), func(x_i))

[A_k, TB, x_i] = fourier_series_expand(func, -1, 1, 11);
assert_equals(trig_eval(A_k, TB, x_i), func(x_i))

[A_k, TB, x_i] = fourier_series_expand(func, -2, 2, 11);
assert_equals(trig_eval(A_k, TB, x_i), func(x_i))

[A_k, TB, x_i] = fourier_series_expand(func, 0, 2, 11);
assert_equals(trig_eval(A_k, TB, x_i), func(x_i))

[A_k, TB, x_i] = fourier_series_expand(func, 0.01, 2.5, 141);
assert_equals(trig_eval(A_k, TB, x_i), func(x_i))
%plot(x_i, func(x_i), x_i, trig_eval(TB, x_i))



func = @(x)(exp(-abs(x)));
[A_k, TB, x_i] = fourier_series_expand(func, -1, 1, 20);
assert_equals(trig_eval(A_k, TB, x_i), func(x_i))

[A_k, TB, x_i] = fourier_series_expand(func, -1, 1, 20, 'symmetry', 'even');
assert_equals(trig_eval(A_k, TB, x_i), func(x_i))


func = @(x)(sign(x).*(exp(-abs(x))-1));
[A_k, TB, x_i] = fourier_series_expand(func, -1, 1, 100, 'symmetry', 'odd');
x_i = x_i(2:end);
assert_equals(trig_eval(A_k, TB, x_i), func(x_i))
