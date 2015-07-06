function [TB, r_k]=kd_fourier(func, pos, varargin)
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
        [S_k, TP] = power_spectrum_by_density(func, L, K, d);
    else
        [S_k, TP] = power_spectrum_by_fft(func, L, K, d);
    end
    if autoenlarge && sum(abs(S_k))-1>1e-7
        L = L * 2;
        strvarexpand('Enlarging to $L(1)$');
    else
        break
    end
end
S_k = S_k';

w_k = TP{1};
p_k = TP{2};
[S_k, w_k, p_k] = sort_spectrum(S_k, w_k, p_k);
[S_k, w_k, p_k] = limit_spectrum(S_k, w_k, p_k, ratio);
assert(all(S_k)>0);

s_k=sqrt(S_k);
% Limit the number of functions
K0 = floor((K-1)/2^d)+1;
if length(s_k)>K0
    s_k = s_k(1:K0);
    w_k = w_k(1:K0,:);
    p_k = p_k(1:K0,:);
end

[s_k, w_k, p_k] = add_sines(s_k, w_k, p_k, d);
TB = {w_k, p_k, s_k};

if ~isempty(pos)
    r_k = trig_basis_eval(TB, pos);
end


function [S_k, w_k, p_k] = sort_spectrum(S_k, w_k, p_k)
% SORT_SPECTRUM Sort the spectral density in descending order
multiplicity = 2.^sum(w_k~=0,2);
Sn_k = S_k.^2 ./ multiplicity;
[~, ind]=sort(Sn_k(2:end), 'descend');
ind = [1; ind+1];

S_k = S_k(ind);
w_k = w_k(ind,:);
p_k = p_k(ind,:);


function [S_k, w_k, p_k] = limit_spectrum(S_k, w_k, p_k, ratio)
% LIMIT_SPECTRUM Limit the number of functions such that RATIO percent of
% the spectrum is covered.
r = cumsum(S_k);
max_ind = find(r>=ratio, 1, 'first');
if isempty(max_ind)
    warning('sglib:kd_fourier', 'Not enough basis functions for given ratio. Increase max_funcs option');
    max_ind = sum(S_k>0);
end
S_k = S_k(1:max_ind);
w_k = w_k(1:max_ind,:);
p_k = p_k(1:max_ind,:);

function [s_k, w_k, p_k] = add_sines(s_k, w_k,p_k, d)
for i=1:d
    ind=w_k(:,i)~=0;

    p_k_s = p_k;
    p_k_s(:,i)=0;

    s_k = [s_k; s_k(ind)];     %#ok<AGROW>
    w_k = [w_k; w_k(ind,:)];   %#ok<AGROW>
    p_k = [p_k; p_k_s(ind,:)]; %#ok<AGROW>
end

function [S_k, TP] = power_spectrum_by_fft(func, L, K, d)
M=2*K+1;
if length(L)==1 && d>1
    L = L * ones(d,1);
end
assert(d==1, 'd>1 does not work yet');
[S_k, TP] = fourier_series_expand(func, -L, L, M, 'symmetry', 'even');


function [S_k, TP] = power_spectrum_by_density(func, L, K, d)
% Hypersphere: V = pi^(d/2) r^d / gamma(d/2+1)
% Radius:      r = (V * gamma(d/2+1))^(1/d) / sqrt(pi)
% K_i = r, V = K * 2^d

%K_i = ceil( 2 * (K / nball_volume(d,1)) ^ (1/d) - 1e-10);
%assert( nball_volume(d, K_i) / 2^d>=K && nball_volume(d, K_i-1) / 2^d<=K)

K_i = ceil( (K / nball_volume(d,1)) ^ (1/d) - 1e-10);
assert( nball_volume(d, K_i-1)<=K && K<=nball_volume(d, K_i)+1 )

I=multiindex(d,K_i,'full', true);
ind = sum(I.^2,2)<=K_i^2;
I = I(ind,:);
w0 = 1./(2*L);
w = binfun(@times, I', w0);
assert(size(I,1)<10000);

S_k = funcall(func, w, d)' * prod(w0);

multiplicity = 2.^sum(w~=0,1)';
S_k = multiplicity .* S_k;
%S_k(2:end) = 2^d * S_k(2:end);
%wp_k = [2*pi*w', repmat(pi/2, size(w'))];
%wp_k = [wp_k(:, 1:2:end), wp_k(:,2:2:end)];
S_k = S_k';
TP = {w', repmat(1/4, size(w'))};
