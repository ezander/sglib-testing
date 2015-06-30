function [A_k, wp_k, x_i] = fourier_series_expand(func, a, b, N, varargin)
% FOURIER_SERIES_EXPAND Short description of fourier_series_expand.
%   [A_K, WP_K, X_I] = FOURIER_SERIES_EXPAND(FUNC, A, B, N, VARARGIN) Long description of fourier_series_expand.
%
% Example (<a href="matlab:run_example fourier_series_expand">run</a>)
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
[symmetry, options]=get_option(options, 'symmetry', 'none');
check_unsupported_options(options);

switch(symmetry)
    case 'even'
        even=true; odd=false;
    case 'odd'
        even=false; odd=true;
    otherwise
        even=true; odd=true;
end

l=b-a;
x_i=linspace(a,b,N+1);
x_i=x_i(1:end-1);
c=fft(funcall(func, x_i));

c=c.*exp(-2*pi*i*a/l*(0:N-1));
a=0;

ar = real(c)/N; 
ai = imag(c)/N;

A_k=[];
w_k=[];
p_k=[];

m=ceil(N/2);
omega = 2*pi/l*(0:m)';
omega_a = a*omega;

% The cosine terms
if even
    if mod(N,2)==0
        A_k=[A_k,  ar(1), 2*ar(2:m), ar(m+1)];
        w_k=[w_k;  omega];
        p_k=[p_k; -omega_a+pi/2];
    else
        A_k=[A_k,  ar(1), 2*ar(2:m)];
        w_k=[w_k;  omega(1:end-1)];
        p_k=[p_k; -omega_a(1:end-1)+pi/2];
    end
end

% The sine terms
if odd
    A_k=[A_k, -2*ai(2:m)];
    w_k=[w_k; omega(2:end-1)];
    p_k=[p_k; -omega_a(2:end-1)];
end

wp_k = [w_k,p_k];
