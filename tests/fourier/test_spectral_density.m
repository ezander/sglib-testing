function test_spectral_density(varargin)
% TEST_SPECTRAL_DENSITY Test whether the spectral density functions integrate to one.

%   Elmar Zander
%   Copyright 2015, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

format compact
format short g

l_c = 0.3;
for i=1:4
    switch i
        case 1
            pow_func = @(w,d)(gaussian_spectral_density(w,l_c,d));
            s = 'gaussian';
        case 2
            pow_func = @(w,d)(exponential_spectral_density(w,l_c,d));
            s = 'exponential';
        case 3
            nu = 1.5;
            pow_func = @(w,d)(matern_spectral_density(nu, w, l_c, d));
            s = 'matern 1.5';
        case 4
            nu = 0.1;
            pow_func = @(w,d)(matern_spectral_density(nu, w, l_c, d));
            s = 'matern 0.1';
    end

    for d=1:8
        I=check_spectral_density(pow_func, d);
        strvarexpand('Cov: $s$, dim=$d$, int=$I$');
    end
end


function I=check_spectral_density(pow_func, d)
integrand = @(t)( pow_func(t,d) .* nball_surface(d, t));
I=integral(integrand, 0, inf);
