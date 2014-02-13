function f=lagrange_newton(xn, fn, x)
% LAGRANGE_NEWTON Short description of lagrange_newton.
%   LAGRANGE_NEWTON Long description of lagrange_newton.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example lagrange_newton">run</a>)
%
% See also

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

F=newton_tableau(xn, fn);
cn = F(1,:);

n=length(fn);

z = ones(size(x));
f = zeros(size(x));
for i=1:n
    f = f + cn(i) * z;
    z = z .* (x-xn(i));
end




