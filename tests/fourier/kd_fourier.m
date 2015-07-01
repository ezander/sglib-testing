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
K0 = floor((K-1)/2^d)+1;
if length(s_k)>K0
    s_k = s_k(1:K0);
    wp_k = wp_k(1:K0,:);
end
[sigma_k, wp_k] = add_sines(s_k, wp_k, d);

if ~isempty(pos)
    r_k = trig_basis_eval(sigma_k, wp_k, pos);
end

function [s_k, wp_k] = sort_spectrum(s_k, wp_k)
[s_k(2:end), ind] = sort(s_k(2:end), 'descend');
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
    wp_k_s = wp_k(2:end,:);
    wp_k_s(:,2*i)=0;
    ind=wp_k_s(:,2*i-1)~=0;

    s_k = [s_k; s_k([false; ind])];
    wp_k = [wp_k; wp_k_s(ind,:)];
end


function [S_k, wp_k] = power_spectrum_by_fft(func, L, K, d)
M=2*K+1;
if length(L)==1 && d>1
    L = L * ones(d,1);
end
assert(d==1, 'd>1 does not work yet');
[S_k, wp_k] = fourier_series_expand(func, -L, L, M, 'symmetry', 'even');

function [S_k, wp_k] = power_spectrum_by_fft(func, L, K, d)
M=2*K+1;
if length(L)==1 && d>1
    L = L * ones(d,1);
end
assert(d==1, 'd>1 does not work yet');
[S_k, wp_k] = fourier_series_expand(func, -L, L, M, 'symmetry', 'even');


function [S_k, wp_k] = power_spectrum_by_density(func, L, K, d)
% Hypersphere: V = pi^(d/2) r^d / gamma(d/2+1)
% Radius:      r = (V * gamma(d/2+1))^(1/d) / sqrt(pi)
% K_i = r, V = K * 2^d
K_i = ceil( 2 * (K * gamma(d/2+1)) ^ (1/d) / sqrt(pi) - 1e-10);
I=multiindex(d,K_i,'full', true);
ind = sum(I.^2,2)<=K_i^2;
I = I(ind,:);
w0 = 1./(2*L);
w = binfun(@times, I', w0);
assert(size(I,1)<10000);

S_k = funcall(func, w, d)' * prod(w0);

if d>=1
    r2 = sum(w.^2, 1);
    S_k = S_k .* (nball_surface(d, sqrt(r2))/2^d)';
    %( (2*pi)^(n/2) / gamma(n/2) * r2 .^ (n/2));
end


S_k(2:end) = 2^d * S_k(2:end);
wp_k = [2*pi*w', repmat(pi/2, size(w'))];
wp_k = [wp_k(:, 1:2:end), wp_k(:,2:2:end)];

