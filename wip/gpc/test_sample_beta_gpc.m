%%
clear
dist_func={@beta_stdnor, {4, 2}, {2, 3}};
[a_i_alpha, I] = pce_expand_1d(dist_func, 4)
V = {'H', I}

N=100000;
xi=gpcgerm_sample(V, N);
y=gpc_evaluate(a_i_alpha, V, xi);
kernel_density(y)


%%
clear

dist_name = 'beta';
params = {4, 2};
dist = gendist_create(dist_name, params);
dist = gendist_fix_moments(dist, 3.2, 0.24);
x=linspace(-1,5);
plot(x,gendist_pdf(x, dist))
hold all
dist_func={@gendist_stdnor, {dist}, {2}}

for p=1:6
    [a_i_alpha, I] = pce_expand_1d(dist_func, p)
    V = {'H', I}

    N=100000;
    xi=gpcgerm_sample(V, N);
    y=gpc_evaluate(a_i_alpha, V, xi);
    empirical_density(y);
    legend('pdf', '1', '2', '3', '4', '5', '6' )
end
