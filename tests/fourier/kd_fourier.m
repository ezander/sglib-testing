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

% Determine dimension and length scale
d = size(pos,1);
if size(pos,2)==1
    L = pos;
    pos = [];
else
    L = max(pos, [], 2) - min(pos, [], 2);
end
if length(L)==1 && d>1
    L = repmat(L, d, 1);
end
    
% Extend the length scale until the spectrum is "ok" (i.e. for the FFT
% based method that all the S_k are non-negative and for the spectral
% density based method that the S_k's nearly add up to one).
% TODO: some less heuristical approach is needed here
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

% Extract the spectrum, sort it (according to amplitude times
% multiplicity), and limit the number of basis functions
w_k = TP{1};
p_k = TP{2};
[S_k, w_k, p_k] = sort_spectrum(S_k, w_k, p_k);
[S_k, w_k, p_k] = limit_spectrum(S_k, w_k, p_k, ratio, K);
assert(all(S_k)>0);

% For the field itself we have to take the square root and add the sin
% terms
s_k=sqrt(S_k);
[s_k, w_k, p_k] = add_sines(s_k, w_k, p_k, d);
TB = {w_k, p_k, s_k};

% Evaluate at given positions, if necessary
if nargout>=2
    r_k = trig_basis_eval(TB, pos);
end


function [S_k, w_k, p_k] = sort_spectrum(S_k, w_k, p_k)
% SORT_SPECTRUM Sort the spectral density in descending order
multiplicity = 2.^sum(w_k~=0,2);
Sn_k = S_k ./ multiplicity;
[~, ind]=sort(Sn_k(2:end), 'descend');
ind = [1; ind+1];

S_k = S_k(ind);
w_k = w_k(ind,:);
p_k = p_k(ind,:);


function [S_k, w_k, p_k] = limit_spectrum(S_k, w_k, p_k, ratio, K)
% LIMIT_SPECTRUM Limit the number of functions such that RATIO percent of
% the spectrum is covered and to max K functions.
r = cumsum(S_k);
max_ind = find(r>=ratio, 1, 'first');
if isempty(max_ind)
    warning('sglib:kd_fourier', ...
            ['Not enough basis functions for given variance ratio (%g). \n '...
            '%d basis functions give only a ratio of %g. '...
            'Increase max_funcs option.'], ...
            ratio, K, max(r));
    max_ind = sum(S_k>0);
end
S_k = S_k(1:max_ind);
w_k = w_k(1:max_ind,:);
p_k = p_k(1:max_ind,:);

% Limit the number of functions by K
multiplicity = 2.^sum(w_k~=0,2);
max_ind = min(max_ind, sum(cumsum(multiplicity)<=K));
if length(S_k)>max_ind
    S_k = S_k(1:max_ind);
    w_k = w_k(1:max_ind,:);
    p_k = p_k(1:max_ind,:);
end


function [s_k, w_k, p_k] = add_sines(s_k, w_k,p_k, d)
% ADD_SINES Add sin basis functions 
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
% POWER_SPECTRUM_BY_DENSITY Compute the power spectrum from the spectral
% density.

% Compute the radius of an n-ball containing K functions.
K_i = ceil( (K / nball_volume(d,1)) ^ (1/d) - 1e-10);
assert( nball_volume(d, K_i-1)<=K && K<=nball_volume(d, K_i)+1 )

% Compute all multindices of frequencies in this n-ball
I=multiindex(d,K_i,'full', true);
ind = sum(I.^2,2)<=K_i^2;
I = I(ind,:);
w0 = 1./(2*L);
w = binfun(@times, I', w0);

% Compute the sprectral density at those frequencies
S_k = funcall(func, w, d)' * prod(w0);

% Adjust for multiplicities (because terms with non-zero frequencies are
% accounted for only once, but really appear two-times or four-times or
% eight-times in the power spectrum...)
multiplicity = 2.^sum(w~=0,1)';
S_k = multiplicity .* S_k;

% Adjust to fit the usual data structures for trig representations
S_k = S_k';
TP = {w', repmat(1/4, size(w'))};
