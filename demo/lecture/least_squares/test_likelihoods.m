function test_likelihoods

%rand_seed(10);

N = 20;
data1 = randn(N, 1)+0.5;
data2 = data1; data2(1)=30;

mu = 0.5;
sigma = 1; % known

mu = linspace(0.01, 2)
data = data2;
for i=1:100
    L1(i) = likelihood_normal(data, mu(i), sigma);
    L2(i) = likelihood_laplace(data, mu(i), sigma);
end
clf

plot(mu,L1/max(L1), mu, L2/max(L2))
grid on


function logL=log_likelihood_normal(m, s)
logL = log(normal_pdf(data, mu, sigma));
logL = sum(logL);


function L = likelihood_normal(data, mu, sigma)
L = normal_pdf(data, mu, sigma);
L = prod(L);

function L = likelihood_laplace(data, mu, sigma)
L = exponential_pdf(abs(data-mu), 1/sqrt(2));
L = prod(L);







