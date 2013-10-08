function [x,w]=newton_cotes_rule(n, open)
% NEWTON_COTES_RULE Compute points and weights of the trapezoidal rule.
%   [X,W]=NEWTON_COTES_RULE(N) computes the points and weights of the
%   closed Newton-Cotes rule with N+1 integration points returning the
%   points in the 1 x N+1 array X and the weights in the N+1 x 1 array w.
%
%   [X,W]=NEWTON_COTES_RULE(N, TRUE) computes the points and weights of the
%   open Newton-Cotes rule with N-1 integration points returning the
%   points in the 1 x N-1 array X and the weights in the N-1 x 1 array w.
%
% References
%   [1] Germund Dalquist, Åke Björck: "Numerical Methods in Scientific
%       Computing, Volume 1", SIAM,
%       http://www.siam.org/books/ot103/OT103%20Dahlquist%20Chapter%205.pdf 
%   [2] https://en.wikipedia.org/wiki/Newton-Cotes_formulas
%
% Example (<a href="matlab:run_example newton_cotes_rule">run</a>)
%   [x,w] = newton_cotes_rule(20);
%   fprintf('Integral of cosine from -1 to 1:\n');
%   fprintf('I_ex     = %g\n', sin(1) - sin(-1))
%   fprintf('I_nc5    = %g\n', cos(x) * w)
%   fprintf('I_matlab = %g\n', quad(@cos, -1, 1))
%
% See also INTEGRATE_1D, GAUSS_LEGENDRE_RULE, TRAPEZOIDAL_RULE



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

if nargin<2
    open = false;
end

if ~open
    % xi = a + i (b − a) / n, i=0..n
    x = linspace(-1, 1, n+1);
else
    % xi = a + i (b − a) / n, i=1..n-1
    x = linspace(-1, 1, n+1);
    x = x(2:end-1);
end
%w = nc_weights_vander(x);
w = nc_weights_lagrange(x);

function w=nc_weights_vander(x)
%% NC_WEIGHTS_VANDER Compute weights via Vandermonde matrix
% int_-1^1 x^n = (1^(n+1) - (-1)^(n+1))/(n+1) = sum_i w_i x_i^n

n1=(1:length(x))';
m=(1 - (-1).^n1)./n1;
V=binfun(@power, x, n1-1);
w=V\m;
return

function w=nc_weights_lagrange(x)
%% NC_WEIGHTS_LAGRANGE Compute weights via Lagrange polynomials
w = zeros(size(x'));
n = length(w);
for k=1:ceil(n/2)
    xk = [x(1:k-1) x(k+1:end)];
    p=poly(xk);
    l=p/polyval(p,x(k));
    l_int=polyint(l);
    w([k,n+1-k])=polyval(l_int,1)-polyval(l_int,-1);
end
