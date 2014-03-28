function y=do_measurement(funcs, u)
% DO_MEASUREMENT Perform a measurement using the measurement operator.
%   Y=DO_MEASUREMENT(FUNCS, U) performs a measurement on U given the
%   functions in the cell array FUNCS. This cell array represents the
%   measurement operator M, that will give actual values that "would be
%   measured" if that U represented the true state of the physical system
%   being modelled.
%
%   U represents a collection of states of the system such that each column
%   vector in U is one distinct state. The FUNCS must therefore return a
%   row vector with one value per distinct state. If there are D functions
%   in FUNCS and U is an N x M matrix, the returned value will be a matrix
%   of size D x M.
%
% Example (<a href="matlab:run_example do_measurement">run</a>)
%   funcs = {@(u)(mean(u,1)), @(u)(u(2,:))};
%   u = rand(10,15);
%   y=do_measurement(funcs, u)
%
% See also

%   Elmar Zander
%   Copyright 2014, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

m = length(funcs);
n = size(u,2);
y = zeros(m,n);
for i=1:m
    y(i,:) = funcall(funcs{i}, u);
end
