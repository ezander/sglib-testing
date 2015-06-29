function [sigma_k, wp_k]=kd_fourier(func, pos, varargin)
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

options=varargin2options(varargin, mfilename);
[K,options]=get_option(options, 'max_funcs', 50);
[is_spectral,options]=get_option(options, 'is_spectral', false);
[autoenlarge,options]=get_option(options, 'autoenlarge', true);
check_unsupported_options(options);

d = size(pos,1);
if size(pos,2)==1
    L = pos;
    pos = [];
else
    L = max(pos, [], 2) - min(pos, [], 2);
end
    

while true
    if is_spectral
        [S_k, wp_k] = power_spectrum_by_density(func, L, K);
    else
        [S_k, wp_k] = power_spectrum_by_fft(func, L, K);
    end
    if autoenlarge && sum(S_k)-1>1e-7
        L(1) = L(1) * 1.5;
        strvarexpand('Enlarging to $L(1)$');
    else
        break
    end
end

s_k=real(sqrt(S_k));
s_k=reshape(s_k,[],1);
sigma_k=[s_k(1); s_k(2:end); s_k(2:end)];
wp_k = [wp_k; wp_k(2:end,1), zeros(K,1)];

if ~isempty(pos)
    r_k = sin_basis_eval(wp_k, pos);
    u_k = binfun(@times, sigma_k, r_k);
end

function [S_k, wp_k] = power_spectrum_by_fft(func, L, K)
M=2*K+1;
[S_k, wp_k] = myfft(func, -L(1), L(1), M, 'symmetry', 'even');


function [S_k, wp_k] = power_spectrum_by_density(func, L, K)
w = (0:K)'/(2*L);
S_k = funcall(func, w) / (2*L);
S_k(2:end) = 2*S_k(2:end);
wp_k = [2*pi*w, repmat(pi/2, size(w))];
