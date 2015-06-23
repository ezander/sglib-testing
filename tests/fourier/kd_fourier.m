function [sigma_k, wp_k]=kd_fourier(func, I, K, varargin)
% KD_FOURIER Short description of kd_fourier.
%   KD_FOURIER(VARARGIN) Long description of kd_fourier.
%
% Example (<a href="matlab:run_example kd_fourier">run</a>)
%
% See also

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

M=2*K+1;
L=I(:,2)-I(:,1);
while true
    [S_k, wp_k] = myfft(func, -L(1), L(1), M, 'symmetry', 'even');
    if norm(S_k(S_k<0)) > 1e-10*norm(S_k)
        %break;
        L(1) = L(1) * 1.5;
        strvarexpand('Enlarging to $L(1)$');
    else
        break;
    end
end

s_k=reshape(sqrt(S_k),[],1);
s_k=real(s_k);
sigma_k=[s_k(1); s_k(2:end); s_k(2:end)];
wp_k = [wp_k; wp_k(2:end,1), zeros(K,1)];
