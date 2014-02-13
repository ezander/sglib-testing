function F=newton_tableau(x, f)
% NEWTON_TABLEAU Short description of newton_tableau.
%   NEWTON_TABLEAU Long description of newton_tableau.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example newton_tableau">run</a>)
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

f=f(:);
x=x(:);

n=size(f,1);
F=zeros(n,n);

F(:,1)=f;
for i=1:n-1
    F(1:n-i, i+1) = (F(1:n-i,i) - F(2:n-i+1,i)) ./ (x(1:n-i) - x(i+1:n));
end
