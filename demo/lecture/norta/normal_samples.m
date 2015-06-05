function xi = normal_samples(d, N, qmc)
% NORMAL_SAMPLES Generate samples from normal distribution

if nargin<3
    qmc = false;
end

if qmc
    xi = normal_invcdf(halton_sequence(N, d)');
else
    xi = randn(d, N);
end
