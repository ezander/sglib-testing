function [A_k, wp_k, x_i] = myfft(func, a, b, N, varargin)
% MYFFT Short description of myfft.
%   MYFFT(VARARGIN) Long description of myfft.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example myfft">run</a>)
%
% See also

%   Elmar Zander
%   Copyright 2015, Inst. of Scientific Computing, TU Braunschweig
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
c=fft(func(x_i));
ar = real(c); ai = imag(c);

A_k=[];
w_k=[];
p_k=[];

% The cosine terms
m=ceil(N/2);
if even
    if mod(N,2)==0
        A_k=[A_k; [ar(1) 2*ar(2:m) ar(m+1)]/N];
        w_k=[w_k; 2*pi*(0:m)'/l];
        p_k=[p_k; -2*pi*(0:m)'*a/l+pi/2];
    else
        A_k=[A_k; [ar(1) 2*ar(2:m)]/N];
        w_k=[w_k; 2*pi*(0:m-1)'/l];
        p_k=[p_k; -2*pi*(0:m-1)'*a/l+pi/2];
    end
end

% The sine terms
if odd
    A_k=[A_k, -2*ai(2:m)/N];
    w_k=[w_k; 2*pi*(1:m-1)'/l];
    p_k=[p_k; -2*pi*(1:m-1)'*a/l];
end

wp_k = [w_k,p_k];
