function x0=fzeros_mult( func, x0, x1, varargin)
% FZEROS_MULT Tries to find all roots of a function in some interval.
%   X0=FZEROS_MULT( FUNC, X0, X1, VARARGIN) tries to find all roots in the
%   interval [X0, X1]. 
%
% Options
%   N: {1000}
%      Initial numbers evaluation points in the interval [x0, x1]
%   tol: {1e-10}
%      Tolerance for call to FZERO.
%
% Note 1: The function may have poles in the interval. Those are ignored
%   for determination of the zeros (see Example).
%
% Note 2: This method is not efficient if function evaluation is expensive.
%   That was not the point of this function. In the first pass it scans the
%   function on a pretty fine mesh, making lots of function evaluations.
%
% Example (<a href="matlab:run_example fzeros_mult">run</a>)
%   % The function below has poles at (n+1/2)*pi, which are not detected
%   % by fzeros_mult
%   fun1=@(x)(1-x.*tan(x));
%   r=fzeros_mult( fun1, 0, 7*pi, 'N', 100 );
%   x=linspace(0, 7*pi);
%   plot(x, fun1(x), 'b-', r, fun1(r), 'rx'); grid on;
%
% See also FZERO

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
[N,options]=get_option(options,'N',1000);
[tol,options]=get_option(options,'tol',1e-10);
check_unsupported_options(options,mfilename);


% evaluate at N points between x0 and x1
x=linspace(x0,x1,N);
y=funcall(func,x);

% find local minima of abs(y), the approach is better than checking sign
% changes of y, if there are discontinuities in the function.
v=abs(y);
s=sign(v(2:end)-v(1:end-1));
k=find(s(2:end)==1 & s(1:end-1)==-1) + 1;

% now from those local minima determin whether there is a sign change and
% keep only those
i1 = find(sign(y(k-1))~=sign(y(k)));
i2 = find(sign(y(k))~=sign(y(k+1)));
k1 = [k(i1)-1, k(i2)];
k2 = [k(i1), k(i2)+1];

x0 = zeros(1,length(k1));
for i=1:length(k1)
    opts=optimset( 'display', 'off', 'TolFun', tol );
    [x0(i),fval,flag]=fzero( @funcall_funfun, x([k1(i), k2(i)]), opts, func );
end
x0=sort(x0);
