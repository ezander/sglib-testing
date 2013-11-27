function x=create_rand_set(N, p, f)

x = randn(N,1);

ind = (rand(N,1)<p);
y = 1+(f-1)*rand(N,1);
x(ind) = x(ind) .* y(ind);
