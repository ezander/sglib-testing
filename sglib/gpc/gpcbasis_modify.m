function V_un = gpcbasis_modify(V_u, varargin)
% GPCBASIS_MODIFY Modifies a GPC basis.
%   V_UN = GPCBASIS_MODIFY(V_U, VARARGIN) Long description of gpcbasis_modify.
%
% Example (<a href="matlab:run_example gpcbasis_modify">run</a>)
%
% See also

%   Elmar Zander
%   Copyright 2015, Insititue of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

V_un = gpcbasis_create(V_u{1}, varargin{:});
