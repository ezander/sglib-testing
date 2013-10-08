function w=interpolatory_weights(x)
%%


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

w = zeros(size(x'));
for k=1:length(w)
    xk = x;
    xk(k)=[];
    p=poly(xk);
    q=polyint(p/polyval(p,x(k)));
    w(k)=polyval(q,1)-polyval(q,-1);
end
