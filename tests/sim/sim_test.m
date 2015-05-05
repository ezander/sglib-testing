% How to use this stuff

kappa = SimParameter('kappa', LogNormalDistribution(1, 2));
mass = SimParameter('mass', NormalDistribution(1, 2));
damping = SimParameter('damping', ExponentialDistribution(32));
time = SimParameter('time', UniformDistribution(3,5));

% sim = SimInterface(@spring_func);
% sim.add_parameter(kappa);
% sim.add_parameter(mass);
% sim.add_parameter(damping);
% sim.add_parameter(time);


param_set = SimParamSet();
param_set.add_parameter(kappa);
param_set.add_parameter(mass);
param_set.add_parameter(damping);
param_set.add_parameter(time);


time.set_to_mean();
mass.set_fixed(5);


xi = param_set.generate_samples(3, 'mode', 'mc');
% should give something like:
% [ 
%     0.234234, 5, 0.56464, 4;
%     0.345344, 5, 4.34534, 4;
%     1.123123, 5, 49058.4, 4;]
% maybe transposed??? I think so



for i=1:3
    spring_func(xi(:,i));
end

xi = param_set.generate_integration_points(5, 'grid', 'smolyak');
