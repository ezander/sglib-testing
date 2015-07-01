function V=nball_volume(n, r)
% NBALL_VOLUME Short description of nball_volume.
%   NBALL_VOLUME(VARARGIN) Long description of nball_volume.
%
% Example (<a href="matlab:run_example nball_volume">run</a>)
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

if nargin<2
    r = 1;
end

V = pi.^(n/2) ./ gamma(n/2+1) .* (r.^n);
