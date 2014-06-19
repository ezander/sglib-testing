function [state, info] = undamped_spring_init( varargin )
% ELECTRICAL_NETWORK_INIT Initialises the structure that keeps the internal
% state of the electrical network example.
%
% Example:
%    state = electrical_network_init('R', 200, 'f0', 5);

options = varargin2options(varargin);
[x0, options] = get_option(options, 'x0', 1);
[v0, options] = get_option(options, 'v0', 0);
[T, options] = get_option(options, 'T', 10);
check_unsupported_options(options, mfilename);

% store everything in the state variable
info = struct();
info.num_params = 2; % m, k
info.num_vars = 2;   % x, v

if nargout==2
    state = struct();
else
    state = info;
end
state.u0 = [x0; v0];
state.T = T;
state.d = 0;
