function [y_j_i] = trig_eval(A_j_k, wp_k_l, x_l_i)
% TRIG_EVAL Evaluate the trigonometric expansion of a function.
%   [Y_J_I] = TRIG_EVAL(A_J_K, WP_K_L, X_L_I) evaluate the function
%   represented by the amplitude array A_J_K and the frequency/phase array
%   WP_K_L at the point X_L_I. The output is given by 
%      Y_J_I = \sum_K A_J_K * \sum_L sin(WP_K_L(2K-1,L)*X_L_I +
%      WP_K_L(2K,L))
%
% Example (<a href="matlab:run_example trig_eval">run</a>)
%
% See also TRIG_BASIS_EVAL, FOURIER_SERIES_EXPAND

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

if iscell(wp_k_l)
    TB = wp_k_l;
    y_k_i = trig_basis_eval(TB, x_l_i);
else
    y_k_i = trig_basis_eval([], wp_k_l, x_l_i);
end
    
y_j_i = A_j_k * y_k_i;
