function density_plot(dist_or_samples, varargin)
% DENSITY_PLOT Short description of density_plot.
%   DENSITY_PLOT Long description of density_plot.
%
% Options
%
% References
%
% Notes
%
% Example (<a href="matlab:run_example density_plot">run</a>)
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

options=varargin2options(varargin);
[do_hold, options]=get_option(options,'hold',false);
% unsupported options are checked inside the plotting subfunctions

%check_type(dist_or_samples, {'cell', 'double'}, false, 'dist_or_samples', mfilename);

if do_hold
    was_hold=ishold();
    hold('all');
end

if iscell(dist_or_samples)
    dist_density_plot(dist_or_samples, options);
else
    sample_density_plot(dist_or_samples, options);
end

if do_hold && ~was_hold
    hold('off');
end

function dist_density_plot(dist, options)
[q0,options] = get_option(options, 'q0', 0);
[q1,options] = get_option(options, 'q1', 1);
[dq,options] = get_option(options, 'dq', 0.02);
[N,options] = get_option(options, 'N', 100);
[ext,options] = get_option(options, 'ext', 0.04);
[type,options] = get_option(options, 'type', 'pdf');
check_unsupported_options(options, [mfilename, '/dist'])

x0 = gendist_invcdf(q0, dist);
if ~isfinite(x0)
    x0 = gendist_invcdf(q0 + dq, dist);
end
x1 = gendist_invcdf(q1, dist);
if ~isfinite(x1)
    x1 = gendist_invcdf(q1 - dq, dist);
end
x = [x0, x1];
if ext~=0
    dx = 0.5 * ext * (x1 - x0);
    x0 = x0 - dx;
    x1 = x1 + dx;
end
x = unique([x, linspace(x0, x1, N)]);
switch(type)
    case 'pdf'
        plot(x, gendist_pdf(x, dist));
    case 'cdf'
        plot(x, gendist_cdf(x, dist));
    case {'both', 'pdf/cdf', 'cdf/pdf'}
        plot(x, gendist_pdf(x, dist), x, gendist_cdf(x, dist));
end        


function sample_density_plot(x, options)
[type, options]=get_option(options, 'type', 'hist');
[n, options]=get_option(options, 'n', 100);
[kde_sig, options]=get_option(options, 'kde_sig', []);
[rug, options]=get_option(options, 'rug', false);
[max_rug, options]=get_option(options, 'max_rug', inf);
check_unsupported_options(options, [mfilename, '/sample']);

switch(type)
    case 'hist'
        kernel_density(x, n, kde_sig);
        %error('foo');
    case 'empirical'
        error('foo');
    case 'kernel'
        kernel_density(x, n, kde_sig);
    otherwise
        error('foo');
end
        
if rug
    if length(x)>max_rug
        x = x(1:max_rug);
    end
    rug_plot(x);
end

