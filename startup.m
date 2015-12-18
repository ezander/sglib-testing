function startup
% STARTUP Called automatically by Matlab at startup.
%   STARTUP gets automatically called by Matlab if it was started from THIS
%   directory (i.e. the sglib home directory). STARTUP just delegates the
%   work SGLIB_STARTUP, so that SGLIB_STARTUP can be called from anywhere
%   else without interfering with any other startup script that might be on
%   the path.
%
% Example (<a href="matlab:run_example startup">run</a>)
%   startup
%
% See also

%   Elmar Zander
%   Copyright 2006, Institute of Scientific Computing, TU Braunschweig.
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.


% We do the real startup in a file with a special name (SGLIB_STARTUP) so
% the user can run it individually without any startup on the path
% interfering with this one
basepath=fileparts( mfilename('fullpath') );
run( fullfile( basepath, 'sglib', 'startup' ) );

% need to get the basepath again, since sglib_startup resets it
basepath=fileparts( mfilename('fullpath') );
addpath( fullfile(basepath, 'data') );
addpath( genpath(fullfile(basepath, 'proposed')) );

addpath(fullfile(basepath, 'demo', 'models'));
%addpath(fullfile(basepath, 'demo', 'methods'));
%addpath(fullfile(basepath, 'demo', 'tutorial', 'tools'));

%addpath(fullfile(basepath, 'tests', 'objects'));

