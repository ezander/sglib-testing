function S=matern_spectral_density(nu, xi, l_c, d)
% MATERN_SPECTRAL_DENSITY Short description of matern_spectral_density.
%   MATERN_SPECTRAL_DENSITY(VARARGIN) Long description of matern_spectral_density.
%
% References:
%   [1] C. E. Rasmussen & C. K. I. Williams, Gaussian Processes for Machine
%       Learning, the MIT Press, 2006, ISBN 026218253X.
%       http://www.gaussianprocess.org/gpml/chapters/RW4.pdf‎
%
% Example (<a href="matlab:run_example matern_spectral_density">run</a>)
%
% See also MATERN_COVARIANCE

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

check_boolean(nu>0, 'input argument "nu" must be positive', mfilename);

r2 = sum(xi.^2, 1);

alpha = 2^d * pi^(d/2) * gamma(nu + d/2) * (2*nu)^(nu);
alpha = alpha / (gamma(nu) * l_c^(2*nu));

S = alpha * (2*nu/l_c^2 + 4*pi^2*r2) .^ -(nu + d/2);
