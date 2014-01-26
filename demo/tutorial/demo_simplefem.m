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
% put into `diffusion_1d_solve` which takes two parameters. Here is an
% example solving the problem for a_1=2 and a_2=100.

a1 = 2; a2 = 100;
[u, a, pos]=diffusion_1d_solve(a1, a2);
plot(pos, u); 


%% Stochastic description of the parameters
% Now suppose the parameters a_1 and a_2 are uncertain and described by
% by random variables. Suppose a_1 follows a Beta[1.2,2] distribution,
% which is shift and scaled such that it has support [0.5, 5]. In sglib
% this can be easily specified within the gendist frame work using the
% gendist_create and gendist_fix_bounds functions (see demo_distributions).

dist_a1 = gendist_create('beta', {1.2, 2});
dist_a1 = gendist_fix_bounds(dist_a1, 0.5, 5);
x = linspace(0, 7);
plot(x, gendist_pdf(x, dist_a1), x, gendist_cdf(x, dist_a1));
legend( 'pdf a1', 'cdf a1');

%%
% The second parameter a_2 follows a shifted and scaled Beta[0.3, 0.6]
% distribution with support on [50, 150]. The distribution of a_2 has
% sharp singularities at both ends.

dist_a2 = gendist_create('beta', {0.6, 0.3});
dist_a2 = gendist_fix_bounds(dist_a2, 50, 150);
x = linspace(0, 200, 1000);
plot(x, gendist_pdf(x, dist_a2), x, gendist_cdf(x, dist_a2));
legend( 'pdf a2', 'cdf a2');

%% GPC approximation of the parameters
% For some methods we could directly use the distributions given above but
% here we want to use their gpc approximation.

%%
% For $a_1$ the distribution looks quite a bit like a (Wigner) semicircle
% distribution for which we can use the combination Semicircle/ChebyshevU,
% specified in sglib by the letter 'u' (small 'u' for the normalised
% polynomials). 
% 
p1 = 2;
V1 = gpcbasis_create('u', 'p', p1);
[x,w]=gpc_integrate([], V1, 10);

psi_k_alpha = gpcbasis_evaluate(V1,x, 'dual', true);
fun_k = gendist_invcdf(gpcgerm_cdf(V1, x), dist_a1{:});
a1_alpha = w'*binfun(@times, psi_k_alpha, fun_k');

% 
% Plot a kernel density estimate of the gpc approximation of a1 and compare
% to the true distribution
%(TODO: gpc_sample should be renamed gpcgerm_sample, and gpc_sample could
% take a parameter a_i_alpha, there should also be a method for
% distribution plotting, rug_plot should determine m automatically) 
a1_samples = gpc_evaluate(a1_alpha, V1, gpc_sample(V1, 100000)); 
kernel_density(a1_samples, 100); hold all;
x=linspace(0,7);
plot(x, gendist_pdf(x, dist_a1{:})); hold off;
legend('gpc approx. (kde)', 'exact density');
rug_plot(a1_samples(1:300));



%% 
%For $a_2$ the distribution looks quite a bit like an arcsine
% distribution for which we can use the combination Arcsine/ChebyshevT,
% specified in sglib by the letter 't'.
p2 = 4;
V2 = gpcbasis_create('u', 'p', p2);
[x,w]=gpc_integrate([], V2, 10);

psi_k_alpha = gpcbasis_evaluate(V2,x, 'dual', true);
fun_k = gendist_invcdf(gpcgerm_cdf(V2, x), dist_a2{:});
a2_alpha = w'*binfun(@times, psi_k_alpha, fun_k');

% 
% Plot a kernel density estimate of the gpc approximation of a2 and compare
% to the true distribution
a2_samples = gpc_evaluate(a2_alpha, V2, gpc_sample(V2, 100000));
kernel_density(a2_samples, 100, 0.01); hold all;
x=linspace(0,2.5);
plot(x, gendist_pdf(x, dist_a2{:})); hold off;
legend('gpc approx. (kde)', 'exact density');
rug_plot(a2_samples(1:300), 3);
