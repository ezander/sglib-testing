function unittest_fzero_illinois
% UNITTEST_FZERO_ILLINOIS Test the FZERO_ILLINOIS function.
%
% Example (<a href="matlab:run_example unittest_fzero_illinois">run</a>)
%   unittest_fzero_illinois
%
% See also FZERO_ILLINOIS, MUNIT_RUN_TESTSUITE 

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

munit_set_function( 'fzero_illinois' );

clc
format compact
format short g
func = @(x)(cos(x)-0.8);
xa = 0;
xb = pi/2;
disp(' ')
fzero_illinois(func, xa, xb, 'maxiter', 100)

% disp(' ')
% fzero_illinois(func, xa, xb, 'maxiter', 100, 'illinois', false)
% 
% disp(' ')
% fzero_illinois(func, xa, xb, 'maxiter', 100, 'bisect', true)


func2 = @(x)([cos(x)-0.8; cos(x)-0.7]);
xa = 0;
xb = pi/2;
disp(' ')
[x,flag,info]=fzero_illinois(func2, xa, xb, 'maxiter', 100, 'd', [1,0]);
[x, acos(0.8)]
info.gval
info.fval
[x,flag,info]=fzero_illinois(func2, xa, xb, 'maxiter', 100, 'd', [0,1]);
[x, acos(0.7)]
info.gval
info.fval
info

[x,flag,info]=fzero_ridders(func2, xa, xb, 'maxiter', 100, 'd', [0,1]);
[x, acos(0.7)]
info.gval
info.fval
info

return


options=struct();
options=optimset(options, 'display', 'iter');
%options=optimset(options, 'tolx', 1e-4);
options=optimset(options, 'outputfcn', @outfun);
fzero(func, [xa, xb], options);

function stop = outfun(x, info, state)
stop = false;
%disp([info.iteration, info.fval, x]);
if abs(info.fval)<1e-8
    stop=true;
end
