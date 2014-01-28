%% Solving stochastic PDEs
% This demonstration show how to solve a parametric stochastic PDE using
% the GPC methods in sglib. 'Parametric' here means, that the uncertainties
% in PDE are given by a finite set of parameters, and not by a stochastic
% field described by infinitely many stochastic variable, which has to be
% (stochastically) discretised beforehand.




%% The deterministic problem
% Lets look at the following boundary value problem 
%
% $$-\frac{\partial}{\partial x}\left(a(x) \frac{\partial}{\partial x}\right)=1 \textrm{ on } 
% \mathcal{D}=[0,1]$$
% 
% where a(x) is a coefficient field given by 
%
% $$a(x) = \left\{\!\! \begin{array}{ll} a_1, & x<0.5 \\ a_2, & x\geq0.5 \end{array} \right. $$
% 

%%
% The code for solving this problem on a grid with N=101 nodes has been
% put into |diffusion_1d_solve| which takes two parameters. Here is an
% example solving the problem for a_1=2 and a_2=100.

a1 = 2; a2 = 100;
[u, a, pos]=diffusion_1d_solve(a1, a2);
plot(pos, u); 





%% Stochastic description of the parameters
% Now suppose the parameters a_1 and a_2 are uncertain and described by
% by random variables. Suppose a_1 follows a Beta[1.2,2] distribution,
% which is shift and scaled such that it has support [0.5, 5]. In sglib
% this can be easily specified within the gendist frame work, using the
% |gendist_create| and |gendist_fix_bounds| functions (see
% demo_distributions). 

a1_dist = gendist_create('beta', {1.2, 2});
a1_dist = gendist_fix_bounds(a1_dist, 0.5, 5);

%%
% The probability distribution function (pdf) and the cumulative
% distribution function (cdf) can be computed with |gendist_cdf| and
% |gendist_pdf| function 
x = linspace(0, 7);
plot(x, gendist_pdf(x, a1_dist), x, gendist_cdf(x, a1_dist));
legend( 'pdf a1', 'cdf a1');

%%
% The second parameter a_2 follows a shifted and scaled Beta[0.3, 0.6]
% distribution with support on [50, 150]. The distribution of a_2 has
% sharp singularities at both ends.

a2_dist = gendist_create('beta', {0.6, 0.3});
a2_dist = gendist_fix_bounds(a2_dist, 50, 150);

x = linspace(0, 200, 1000);
plot(x, gendist_pdf(x, a2_dist), x, gendist_cdf(x, a2_dist));
legend( 'pdf a2', 'cdf a2');





%% GPC approximation of the parameters
% For some methods we could directly use the distributions given above but
% here we want to use their gpc approximation.

%%
% For $a_1$ the distribution looks quite a bit like a (Wigner) semicircle
% distribution for which we can use the combination Semicircle/ChebyshevU,
% specified in sglib by the letter 'u' (small 'u' for the normalised
% polynomials). With the option 'varerr' we instruct |gpc_param_expand| to
% choose such a degree of the gpc polynomial, that the error in variance is
% below 0.001. The option 'fixvar' rescale the gpc coefficients such that
% the variance is fully matched (possibly incurring a worse match in higher
% moments).
% 
[a1_alpha, V1, err] = gpc_param_expand(a1_dist, 'u', 'varerr', 0.001, 'fixvar', true);

%%
% Compare the first four moments of the true distribution of a1 and its GPC
% approximation

[mean,var,skew,kurt]=gendist_moments(a1_dist);
fprintf('Moments (true):\nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)
[mean,var,skew,kurt]=gpc_moments(a1_alpha, V1);
fprintf('Moments (gpc): \nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)

%
% Plot a kernel density estimate of the gpc approximation of a1 and compare
% to the true distribution
% (there should be a method for distribution plotting) 
a1_samples = gpc_evaluate(a1_alpha, V1, gpcgerm_sample(V1, 100000)); 
kernel_density(a1_samples, 100); hold all;
x=linspace(0,7);
plot(x, gendist_pdf(x, a1_dist{:})); hold off;
legend('gpc approx. (kde)', 'exact density');
rug_plot(a1_samples(1:300), 'color', [0.7,0,0]);

% It can be seen in the plot 

%% 
%For $a_2$ the distribution looks quite a bit like an arcsine
% distribution for which we can use the combination Arcsine/ChebyshevT,
% specified in sglib by the letter 't'.

[a2_alpha, V2] = gpc_param_expand(a2_dist, 't', 'varerr', 0.01);

%% 
% Plot a kernel density estimate of the gpc approximation of a2 and compare
% to the true distribution
a2_samples = gpc_evaluate(a2_alpha, V2, gpcgerm_sample(V2, 100000));
kernel_density(a2_samples, 100, 0.01); hold all;
x=linspace(0,2.5);
plot(x, gendist_pdf(x, a2_dist{:})); hold off;
legend('gpc approx. (kde)', 'exact density');
rug_plot(a2_samples(1:300));







%%
%

%%
a1_dist = gendist_create('beta', {1.2, 2});
a1_dist = gendist_fix_bounds(a1_dist, 0.5, 5);

format compact
format short g
for p=1:10
    [a1_alpha, V1] = gpc_param_expand(a1_dist, 'u', p);
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
    [a2_alpha, V2] = gpc_param_expand(a2_dist, 't', p);
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
