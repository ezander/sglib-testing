%% Probablity distributions in sglib
% This tutorial discribes how to specify and work with probablity
% distributions in sglib. 

%% Types of distributions
% Currenly, sglib defines only continuous probablity distribtions. From
% those the normal, lognormal, beta, exponential and uniform distribution
% are defined. For each of those distributions there are four functions
% that are needed for most applications: computation of moments, the
% probability density function (pdf), the cumulative distribution function
% (cdf) and the inverse cdf or quantile function. Each of those function
% begins with the name of the distribution and then "_moments", "_pdf",
% "_cdf" and "_invcdf". 
%
% To get the first four moments of the Beta distribution with parameters
% a=2 and b=3, you can write

[mn,vr,sk,kr]=beta_moments(2,3);
fprintf('Beta(2,3): mean=%g, var=%g, skewness=%g, kurtosis excess=%g\n', mn, vr, sk, kr);

%%
% For the normal distribution you specify the parameters as mean (mu) and
% standard deviation (sigma, not sigma^2), so the variance for the normal
% distribution with parameters 2 and 3 will be 9: 

[mn,vr,sk,kr]=normal_moments(2,3);
fprintf('N(2,3): mean=%g, var=%g, skewness=%g, kurtosis excess=%g\n', mn, vr, sk, kr);

%%
% In order to get the plot range for a distribution, you can use the
% quantile function. For e.g. for the normal distribution, if you want to
% plot over the typical 95% confidence interval:

a = 0.05;
xm = normal_invcdf([a/2, 1-a/2],0,1);
fprintf('N(2,3): x_min(%g)=%g, x_max(%g)=%g\n', a/2, xm(1), 1-a/2, xm(2));
x = linspace(xm(1), xm(2));
plot(x, normal_pdf(x,0,1));

%% General distributions 
% 
% In general it is awkward if you want to work with a certain distribution,
% to specify its name and the parameters every time. Furthermore, you may
% want to pass the distribution itself as a parameter to other functions.
% This can be done in sglib with the gendist framework. You create a
% discription of a distribution with the method gendist_create, e.g.:

dist = gendist_create('lognormal', {1, 0.7}, 'shift', 1.3);

%%
% This create a lognormal distribution with parameters mu=1 and sigma=2,
% specified as cell array. The optional parameter shift, shifts the
% distribution to the right, so that it's defined on [1.3,inf] now. To see
% that we can plot the probability density function and the cumulative
% distribution function:

x=linspace(0,10);
plot(x, gendist_pdf(x, dist), x, gendist_cdf(x, dist));
legend('pdf', 'cdf');

%%
% There is another option 'scale' for ...

%%
% also there is fix_moments, fix_bounds, sample, 
% fix the stdnor stuff, (define only when useful)


%% Extending by more distributions
%

%% Using distributions from the statistics toolbox
% sglib has been written such that it does not depend on extra toolboxes
% but only on basic matlab commands and facilities. However, if you have
% the statistics toolbox, you can also 





%%
%

%%
a1_dist = gendist_create('beta', {1.2, 2});
a1_dist = gendist_fix_bounds(a1_dist, 0.5, 5);

format compact
format short g
for p=1:10
    [a1_alpha, V1] = gpc_param_expand(a1_dist, 'u', 'p', p);
    a1_alpha
    [mu,var] = gendist_moments(a1_dist)
    [mu,var] = gpc_moments(a1_alpha, V1)
end


%%
a2_dist = gendist_create('beta', {0.6, 0.3});
a2_dist = gendist_fix_bounds(a2_dist, 50, 150);
format compact
format short g
for p=1:10
    [a2_alpha, V2] = gpc_param_expand(a2_dist, 't', 'p', p);
    a2_alpha
    [mu,var] = gendist_moments(a2_dist)
    [mu,var] = gpc_moments(a2_alpha, V2)
end


%% Deeper parameter expansion studies
% Study the dependence on the polynomial base
a1_dist = gendist_create('beta', {1.2, 2});
a1_dist = gendist_fix_bounds(a1_dist, 0.5, 5);
[a1_alpha, V1] = gpc_param_expand(a1_dist, 'p', 'varerr', 0.001)
[m1,m2,m3,m4]=gpc_moments(a1_alpha, V1); [m1,m2,m3,m4]
[m1,m2,m3,m4]=gendist_moments(a1_dist); [m1,m2,m3,m4]





%% Some more interesting stuff
[x,y]=plot_density(a2_dist, 'N', 1000, 'q0', 0.95);
legend_add('exact');
for p=5:20
    %[a2_alpha_p, V2_p] = gpc_param_expand(a2_dist, 't', 'varerr', 0.002);
    [a2_alpha_p, V2_p, varerr] = gpc_param_expand(a2_dist, 't', 'p', p);
    hold all; plot(x, gpc_pdf_1d(a2_alpha_p, V2_p, x)); hold off;
    legend_add(strvarexpand('p=$p$, $length(a2_alpha_p)$ terms, varerr=$varerr$'));
end
%%
[x,y]=plot_density(a1_dist, 'N', 1000, 'q0', 0);
legend_add('exact');
for varerr=10.^-(1:0.5:3)
    %[a2_alpha_p, V2_p] = gpc_param_expand(a2_dist, 't', 'varerr', 0.002);
    [a1_alpha_p, V1_p, varerr1] = gpc_param_expand(a1_dist, 'p', 'varerr', varerr);
    hold all; plot(x, gpc_pdf_1d(a1_alpha_p, V1_p, x)); hold off;
    p=length(a1_alpha_p)-1;
    legend_add(strvarexpand('p=$p$, $length(a1_alpha_p)$ terms, varerr=$varerr$'));
end
a1_samples = gpc_sample(a1_alpha_p, V1_p, 100000);
plot_density(a1_samples, 'type', 'hist', 'hold', true);
legend_add(strvarexpand('kde'));
