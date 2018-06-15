function [x,y]=create_data(type, N, varargin)

options = varargin2options(varargin);

switch type
    case {1, 'line'}
        m = 2.1;
        n = -0.3;
        sigma = get_option(options, 'sigma', 0.2);
        x = 2 + 3*lhs_uniform(N,1);
        y = m*x + n + sigma * rand(N,1);
    otherwise
        error('Unknown type');
end
