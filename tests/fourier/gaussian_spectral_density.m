function S=gaussian_spectral_density(xi, l_c, d)
% GAUSSIAN_SPECTRAL_DENSITY Computes the spectral density for a Gaussian covariance.
%   GAUSSIAN_SPECTRAL_DENSITY(XI, L_C, D) computes ...
%
% Example (<a href="matlab:run_example gaussian_spectral_density">run</a>)
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
S = (l_c^2*pi)^(d/2) * exp(-pi^2 * r2 * l_c^2);
