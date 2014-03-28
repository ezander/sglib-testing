function [ym_beta, V_ym] = create_measurement_gpc(dists, polysys)
m = length(dists);
ym_beta = zeros(0,1);
V_ym = gpcbasis_create('');
for i=1:m
    [yi_beta, V_yi] = gpc_param_expand(dists{i}, polysys{i});
    [ym_beta, V_ym] = gpc_combine_inputs(ym_beta, V_ym, yi_beta, V_yi);
end



