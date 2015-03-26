function A=chopabs( A, delta )
% CHOPABS Replace numbers close to zero with zero.
%   CHOPABS(A) returns A with all numbers in A smaller then 1e-10 in magnitude
%   replaced by 0. CHOPABS(A,DELTA) returns A with all numbers in A smaller
%   then DELTA in magnitude replaced by 0.
%
% Note: this function was named CHOPABS to distinguish it from the CHOP
%   function of the control toolbox and further to make clear that it does
%   absolute chopping; not relative like the control toolbox function does.
%
% Example (<a href="matlab:run_example chopabs">run</a>)
%   A=[ 1 2 1e-11 1e-3 3];
%   disp(A);
%   A=chopabs( A );
%   disp(A);
%
% See also ROUND, CEIL, FLOOR

%   Elmar Zander
%   Copyright 2006, Institute of Scientific Computing, TU Braunschweig.
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

if nargin<2
    delta=1e-10;
end
A(abs(A)<delta)=0;
