function plot_samples(xi, varargin)

options=varargin2options(varargin);
[color,options]=get_option(options,'color','b');
check_unsupported_options(options,mfilename);



plot(xi(1,:), xi(2,:), 'LineStyle', 'none', 'Marker', '.', ...
    'MarkerSize', 0.01, 'LineWidth', 0.01, ...
    'MarkerEdgeColor', color, 'MarkerFaceColor', color); 
