function [y_k_i] = trig_basis_eval(a_k, wp_k_l, x_l_i)
% TRIG_BASIS_EVAL Evaluates sequence of trigonometric basis functions.
%   [Y_K_I] = TRIG_BASIS_EVAL(WP_K_L, X_L_I)  evaluate the trigonometric
%   basis functions represented by the amplitude array A_K and the
%   frequency/phase array WP_K_L at the point X_L_I. The output is given by
%      Y_K_I = \sum_L A_K(K) sin(WP_K_L(2K-1,L)*X_L_I + WP_K_L(2K,L))
%   If the A_K is empty, it is assumed the amplitudes are 1.
%
% Example (<a href="matlab:run_example trig_basis_eval">run</a>)
%
% See also TRIG_EVAL

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

if nargin>2
    tau = 2*pi;
    TB = {wp_k_l(:, 1:2:end)/tau, wp_k_l(:, 2:2:end)/tau, a_k};
    y_k_i = trig_basis_eval_int(TB, x_l_i);
else
    y_k_i = trig_basis_eval_int(a_k, wp_k_l);
end

function [y_k_i] = trig_basis_eval_int(TB, x_l_i)

w_k_l = TB{1};
p_k_l = TB{2};
if length(TB)>=3
    a_k = TB{3};
else
    a_k = [];
end

if isempty(a_k)
    y_k_i = ones(size(w_k_l, 1), size(x_l_i, 2));
else
    a_k = a_k(:);
    y_k_i = repmat(a_k, 1, size(x_l_i, 2));
end

tau = 2*pi;
for l=1:size(x_l_i,1)
    y_k_i = y_k_i .* sin(tau * binfun(@plus, w_k_l(:,l)*x_l_i(l,:), p_k_l(:,l)));
end
