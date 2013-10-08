function rug_plot(x, m, varargin)
% RUG_PLOT Creates a rug plot of the given data.
%   RUG_PLOT(X, M) plots small strips of length M*SCALE at the data point
%   given by X. SCALE is by default 0.04 and M should be set to the maximum
%   probability density of the data given by X (estimated e.g. by a kernel
%   density estimate).
%
% Options
%   scale: 0.04
%     Multiplication factor for M. Can also be negative, then the "fringes"
%     point to the negative direction.
%
% Example (<a href="matlab:run_example rug_plot">run</a>)
%     x = randn(1, 300);
%     [xd, pd] = kernel_density(x);
%     plot(xd, pd);
%     hold on;
%     rug_plot(x, max(pd(:)));
%     hold off;
%
% See also KERNEL_DENSITY, PLOT, EMPIRICAL_DENSITY

%   Elmar Zander
%   Copyright 2013, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

options=varargin2options(varargin);
[scale, options]=get_option(options, 'scale', 0.04);
check_unsupported_options(options, mfilename);

x=x(:);
z = zeros(size(x));
X=reshape([x, x, nan*x]', [], 1);
Y=reshape([z, z+scale * m, nan*x]', [], 1);
line(X, Y, 'Color', 'r')
