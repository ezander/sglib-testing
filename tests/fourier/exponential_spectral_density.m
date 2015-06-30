function S=exponential_spectral_density(xi, l_c, d)
% EXPONENTIAL_SPECTRAL_DENSITY Computes the spectral density for an exponential covariance.
%   EXPONENTIAL_SPECTRAL_DENSITY(XI, L_C, D) computes ...
%
% Example (<a href="matlab:run_example exponential_spectral_density">run</a>)
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

if nargin<3 || isempty(d)
    d = size(xi,1);
end

r2 = sum(xi.^2, 1);
alpha = (2*l_c)^d * pi^((d-1)/2) * gamma((d+1)/2);

S = alpha * (1 + 4*pi^2*l_c^2*r2) .^ -((d+1)/2);
