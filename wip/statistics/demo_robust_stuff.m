N=1000;
p=0.01;
%p=0;
f=7;
x=create_rand_set(N, p, f);
x = x.^2;
max(x)
kernel_density(x);

trimmed_mean(x,0)
trimmed_mean(x,0.1)

