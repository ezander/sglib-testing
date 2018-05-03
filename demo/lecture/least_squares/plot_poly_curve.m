function plot_poly_curve(q, a, b, varargin)

options=varargin2options(varargin);
[color, options]=get_option(options, 'color', 'k');

x = linspace(a, b, 1000);
y = polyval(q, x);
line(x, y, 'color', color);
