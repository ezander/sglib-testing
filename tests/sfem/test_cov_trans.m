dist=gendist_create('beta', {3,2});
x=linspace(-1,2);
plot(x,gendist_pdf(x, dist))
rho_stdnor_func=funcreate(@gendist_stdnor, @funarg, dist);


[rho_k, I]=pce_expand_1d(rho_stdnor_func,3)
hold all
y=pce_evaluate(rho_k, I, x)
plot(x,y)
hold off

