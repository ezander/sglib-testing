function [sigma_k, wp_k, r_k]=kd_fourier(func, pos, varargin)
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
[K,options]=get_option(options, 'max_funcs', 2000);
[ratio,options]=get_option(options, 'ratio', 0.99);
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
        [S_k, wp_k] = power_spectrum_by_density(func, L, K, d);
    else
        [S_k, wp_k] = power_spectrum_by_fft(func, L, K, d);
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
[s_k, wp_k] = sort_spectrum(s_k, wp_k);
[s_k, wp_k] = limit_spectrum(s_k, wp_k, ratio);
[sigma_k, wp_k] = add_sines(s_k, wp_k, d);

if ~isempty(pos)
    r_k = trig_basis_eval(sigma_k, wp_k, pos);
end

function [s_k, wp_k] = sort_spectrum(s_k, wp_k)
[s_k(2:end), ind] = sort(s_k(2:end));
wp_k(2:end,:) = wp_k(ind+1,:);

function [s_k, wp_k] = limit_spectrum(s_k, wp_k, ratio)
r = cumsum(s_k.^2);
max_ind = find(r>=ratio, 1, 'first');
if isempty(max_ind)
    warning('sglib:kd_fourier', 'Not enough basis functions for given ratio. Increase max_funcs option');
else
    s_k = s_k(1:max_ind);
    wp_k = wp_k(1:max_ind,:);
end

function [s_k, wp_k] = add_sines(s_k, wp_k, d)
for i=1:d
    s_k=[s_k(1); s_k(2:end); s_k(2:end)];
    wp_k_s = wp_k(2:end,:);
    wp_k_s(:,2*i)=0;
    wp_k = [wp_k; wp_k_s];
end


function [S_k, wp_k] = power_spectrum_by_fft(func, L, K, d)
M=2*K+1;
if length(L)==1 && d>1
    L = L * ones(d,1);
end
assert(d==1, 'd>1 does not work yet');
[S_k, wp_k] = fourier_series_expand(func, -L, L, M, 'symmetry', 'even');


function [S_k, wp_k] = power_spectrum_by_density(func, L, K, d)
K_i = ceil(K ^ (1/d));
w_i = {};
for i = 1:d
    w_i{i} = (0:K_i)/(2*L(i));
end
w = tensor_mesh(w_i, w_i);

S_k = funcall(func, w, d)' / prod(2*L);
S_k(2:end) = 2^d * S_k(2:end);
wp_k = [2*pi*w', repmat(pi/2, size(w'))];
wp_k = [wp_k(:, 1:2:end), wp_k(:,2:2:end)];

function tensorize(s_k, wp_k, d)
K_i = ceil(K0 ^ (1/d));
S_k = ones(1);
wp_k = zeros(1,0);
K = 1;
for i = 1:d
    w = (0:K_i)'/(2*L(i));
    S_k_i = funcall(func, w', d)' / (2*L(i));
    S_k_i(2:end) = 2*S_k_i(2:end);
    wp_k_i = [2*pi*w, repmat(pi/2, size(w))];

    I1 = repmat((1:K)', K_i, 1);
    I2 = repmat((1:K_i), 1, K);
    S_k = S_k(I1) .* S_k_i(I2);
    wp_k = [wp_k(I1,:), wp_k_i(I2,:)];
    K = K_i * K;
end

