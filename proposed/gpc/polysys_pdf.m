function p=polysys_pdf(sys, xi)
% POLYSYS_PDF Evaluates the PDF for an associated probability measure.
%   P=POLYSYS_PDF(SYS, XI) evaluates the probability distribution function
%   (PDF) of the probability measure associated with the polynomial system
%   SYS at the points given in XI. E.g. POLYSYS_PDF('H', 3) would evaluate
%   the normal distribution at X=3, since the normal distribution is
%   associated with the Hermite polynomials (H).
%
% Example (<a href="matlab:run_example polysys_pdf">run</a>)
%   x=linspace(-3,3,400);
%   plot(x, polysys_pdf('h',x), x, polysys_pdf('p',x), ...
%     x, polysys_pdf('l',x), x, polysys_pdf('t',x), x, polysys_pdf('u',x));
%   legend('Gauss/Hermite (h)', 'Uniform/Legendre (p)', ...
%          'Exponential/Laguerre (l)', 'Arcsine/ChebyshevT (t)', ...
%          'Semicircle/ChebyshevU (u)');
% See also GPC, GPCGERM_PDF

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

switch upper(sys)
    case 'H'
        p = normal_pdf(xi, 0, 1);
    case 'P'
        p = uniform_pdf(xi, -1, 1);
    case 'T'
        % Arcsine distribution with support [-1,1] (which is the same as as
        % Beta(1/2,1/2) distribution with shifted support.
        p = beta_pdf(0.5*(xi+1), 1/2, 1/2);
    case 'U'
        % Wigner semicircle distribution (which is the same as as
        % Beta(3/2,3/2) distribution shift from [0,1] to [-1,1].
        p = beta_pdf(0.5*(xi+1), 3/2, 3/2);
    case 'L'
        % Exponential distribution
        p = exponential_pdf(xi, 1);
    case 'M'
        error('sglib:gpc:polysys', 'There is no distribution associated with the monomials.');
    otherwise
        error('sglib:gpc:polysys', 'Unknown polynomials system: %s', sys);
end
