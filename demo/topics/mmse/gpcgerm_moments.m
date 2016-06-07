function [mean,var,skew,kurt]=gpcgerm_moments( V_r )
% GPCGERM_MOMENTS Short description of gpcgerm_moments.
%   GPCGERM_MOMENTS(VARARGIN) Long description of gpcgerm_moments.
%
% Example (<a href="matlab:run_example gpcgerm_moments">run</a>)
%
% See also

%   Elmar Zander
%   Copyright 2016, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

n = nargout;
m = gpcbasis_info(V_r, 'germ_size');

moments = num2cell(nan(m, 4));
syschars = V_r{1};
oldsyschar = -1;

for i=1:m
    syschar = syschars(min(i, length(syschars)));
    if syschar~=oldsyschar
        dist = polysys_dist(syschar);
        [moments{i,1:n}] = gendist_moments(dist);
    else
        moments(i,:) = moments(i-1,:);
    end
    oldsyschar = syschar;
end    

moments = cell2mat(moments);
moments = mat2cell(moments, m, [1,1,1,1]);
[mean,var,skew,kurt]=deal(moments{:});
