function [ym_beta, V_ym] = create_measurement_gpc(dists, polysys)
% CREATE_MEASUREMENT_GPC Create a GPC expansion for the measurements.
%   [YM_BETA, V_YM] = CREATE_MEASUREMENT_GPC(DISTS, POLYSYS) creates from
%   cell array DISTS containing the distributions of independent
%   measurements a GPC that represents all those distributions as ouput. In
%   the cell array POLYSYS the polynomial system used for each basis random
%   variable in the GPC germ needs to be specified.
%
% Example (<a href="matlab:run_example create_measurement_gpc">run</a>)
%   dists{1} = gendist_create('normal', {10, 0.5});
%   polysys{1} = 'H';
%   dists{2} = gendist_create('uniform', {3.5, 6.5});
%   polysys{2} = 'P';
%   [ym_beta, V_ym] = create_measurement_gpc(dists, polysys)
%   ym_samples = gpc_sample(ym_beta, V_ym, 10000);
%   plot( ym_samples(1,:), ym_samples(2,:), '.');
%
% See also GENDIST_CREATE, GPC_EVALUATE, GPC_INTEGRATE, GPC

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

m = length(dists);
ym_beta = zeros(0,1);
V_ym = gpcbasis_create('');
for i=1:m
    [yi_beta, V_yi] = gpc_param_expand(dists{i}, polysys{i});
    [ym_beta, V_ym] = gpc_combine_inputs(ym_beta, V_ym, yi_beta, V_yi);
end
