%% Solving stochastic PDEs
% This demonstration shows how to solve a parametric stochastic PDE using
% the GPC methods in sglib. 'Parametric' here means, that the uncertainties
% in the SPDE are given by a finite set of parameters, and not by a
% stochastic field described by infinitely many stochastic variables, which
% have to be (stochastically) discretised beforehand.

%% The deterministic problem
% We consider the following boundary value problem 
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
% put into the function |diffusion_1d_complete_solve|, which takes the two
% parameters a1 and a2 as arguments. Here is an example solving the problem
% for a1=2 and a2=100. 

a1 = 2; a2 = 100;
[u, a, pos]=diffusion_1d_complete_solve(a1, a2);
plot(pos, u); 


%% Stochastic description of the parameters
% Now suppose the parameters a1 and a2 are uncertain and described by
% by random variables. To make things a bit interesting, we take two
% different distributions for a1 and a2, with very different shapes and
% very different ranges.
% 
% Suppose therefore, that a1 follows a Beta(1.2,2) distribution, which is
% shift and scaled such that it has support [0.5, 5]. In sglib this can be
% easily specified with the statistics functions, using the
% |gendist_create| and |gendist_fix_bounds| functions (|gendist| stands
% here for general distribution, because the functions work on
% distributions specified by cell array objects and relay the calls to the
% concrete underlying distributions.).  

a1_dist = gendist_create('beta', {1.2, 2});
a1_dist = gendist_fix_bounds(a1_dist, 0.5, 5);

%%
% The probability distribution function (pdf) and the cumulative
% distribution function (cdf) can be computed with |gendist_cdf| and
% |gendist_pdf| function. Both of which are used in the |plot_density|
% function. As can be seen in the plot the distribution of a1 is unimodal,
% but asymmetrical leaning somewhat to the left. 

plot_density(a1_dist, 'type', 'both')
legend( 'pdf a1', 'cdf a1');

%%
% The second parameter a2 follows a shifted and scaled Beta(0.7, 0.6)
% distribution with support on [50, 150]. The distribution of a2 is also
% asymetrical, but it is bimodal and has sharp singularities at both ends.
% In order to better see the singularities in the plot some extra x-values
% to plot in their vicinity are specified with the 'extra_x' option.

a2_dist = gendist_create('beta', {0.7, 0.6});
a2_dist = gendist_fix_bounds(a2_dist, 50, 150);

plot_density(a2_dist, 'type', 'both', 'extra_x', [50.0001, 149.999])
legend( 'pdf a2', 'cdf a2');

%% GPC approximation of the parameters
% For some stochastic methods we could directly use the distributions given
% above, but in this demonstration we want to use their respective GPC
% approximations.

%%
% When we select a basis for the GPC appoximation, it is usually best to
% select one with a distribution, that resembles the distribution of the
% parameter, so that only few terms in the GPC expansion are needed. The
% distribution of the the parameter a1 looks a bit like a (Wigner)
% semicircle distribution,  which has the shape of the upper half  of a
% circle between -1 and 1, normalised to have area 1. 
% The polynomials corresponding to this distribution are the Chebyshev
% polynomials of the second kind, usually denoted by the letter U. 
%
% To specify this combination (Semicircle/ChebyshevU) as GPC basis in
% sglib, we use the letter 'u'. The small 'u' selects the normalised
% polynomials in contrast to the unnormalised polynomials 'U'. With the
% option 'varerr' we instruct |gpc_param_expand| to choose such a degree of
% the gpc polynomial, that the error in variance is below 0.001. The option
% 'fixvar' rescales the GPC coefficients such that the variance is fully
% matched (possibly incurring a worse match in higher moments).

[a1_alpha, V1, err] = gpc_param_expand(a1_dist, 'u', 'varerr', 0.001, 'fixvar', true);

%%
% To check the quality of the approximation, we compare the first four
% moments of the true distribution of a1 with the moments of its GPC
% approximation:

[mean,var,skew,kurt]=gendist_moments(a1_dist);
fprintf('Moments (a1,true):\nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)
[mean,var,skew,kurt]=gpc_moments(a1_alpha, V1);
fprintf('Moments (a1,gpc): \nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)

%%
% Often, it is also a good idea to compare the distributions visually,
% which can be achieved by the |plot_density| function. This function can
% either plot an exact pdf or approximate a pdf from sampled data, using
% kernel density estimates or histograms. So for a1, we first draw some
% samples from the underlying distribution, then plot a histogram, then a
% kernel density estimate of the GPC approximation of a1 using those
% samples and then compare this to the true distribution.
%
% The 'rug' option in the second call, instructs |plot_density| to include a
% rug plot, so that we can see where the samples truly were (because the
% kernel density estimator is smearing out the distribution, sometimes
% giving a false impression of the range of the data).
a1_samples = gpc_sample(a1_alpha, V1, 100000);

plot_density(a1_samples, 'n', 30);
plot_density(a1_samples, 'type', 'kernel', 'rug', true, 'max_rug', inf, 'hold', true);
plot_density(a1_dist, 'hold', true);
legend('GPC histogram', 'GPC kernel density', 'GPC samples', 'Exact pdf');

%%
% The approximation of the shape of distribution is not very good in the
% maximum of the distribution. If you change the 'varerr' to 0.0001 in the
% call to |gpc_param_expand|, you can see that the approximation becomes
% much better, however at the cost of 7 expansion coefficients instead of
% just 3. 

%% 
% The distribution of the parameter a2 resembles quite much an arcsine
% distribution. This distribution goes with the Chebyshev polynomials of
% the first kind as family of orthogonal polynomials, which are commonly
% denoted by the letter 'T'. So in sglib, to get the normalised Chebyshev T
% polynomials with corresponding Arcsine distribution we use the letter
% 't'. 

[a2_alpha, V2] = gpc_param_expand(a2_dist, 't', 'varerr', 0.002, 'fixvar', true);

%%
% Now again we compare the first four moments of the true distribution of
% a2 with its GPC approximation:
[mean,var,skew,kurt]=gendist_moments(a2_dist);
fprintf('Moments (a2,true):\nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)
[mean,var,skew,kurt]=gpc_moments(a2_alpha, V2);
fprintf('Moments (a2,gpc): \nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)

%%
% The moment approximations looked quite well. Now, we plot a kernel
% density estimate of the GPC approximation of a2 and compare that to the
% true distribution. Note that for the kernel density estimation we need to
% specify the with of the kernels (setting 'kde_sig' to 0.5), because the
% default used here (Silverman's rule of thumb) is not appropriate for this
% kind of bimodal distribution.
a2_samples = gpc_sample(a2_alpha, V2, 100000);

plot_density(a2_samples, 'type', 'kernel', 'kde_sig', 0.5, 'rug', true, 'max_rug', 300);
plot_density(a2_dist, 'hold', true);
legend('GPC kernel density', 'GPC samples', 'Exact pdf');


%% Combining the GPC spaces
% The GPC spaces V1 anv V2, constructed for the stochastic discretisation
% of the input parameters a1 and a2, are currently representations of
% disjoint probablity spaces. To use them in any of the further algorithms
% (except maybe Monte Carlo), we need to combine them into a representation
% of a larger space, i.e. the Cartesian product of those spaces with the
% corresponding product measure). This can be achieved by the
% |gpc_combine_inputs| function (which can be used much more generally than
% the name implies.). This will not only form the product space, but also
% map the coefficients to the right places in the new combined coefficient
% field.

[a_i_alpha, V_a] = gpc_combine_inputs(a1_alpha, V1, a2_alpha, V2);

%%
% The object V_a that represents the GPC space, can usually be treated as
% abstract entity and gets just passed around along the coefficients in
% a_i_alpha. If you are interested: here is what is inside. V_a is a cell
% array, whose first entry contains the germ of the space (in this case
% 'ut') and the second entry contains the multiindices of the orthogonal
% polynomials corresponding to the germ.

underline('germ'); disp(V_a{1});
underline('multiindices', 'newlines', 1); disp(V_a{2});

%%
% So here e.g. the fourth basis function is given by the multiindex (0,3)
% and thus corresponds to the polynomial u_0(r)t_3(s). With the
% |gpcbasis_polynomials| function we can look at the polynomials that form
% the basis.

cellfun(@(p)fprintf('%s\n',p), ...
    gpcbasis_polynomials(V_a, 'symbols', 'r,s'))

%%
% The coefficients in a_i_alpha are a 2 dimensional array, in which the
% first row corresponds to a1 and the second row to a2. The columns in the
% coefficients are related to the correspondings rows in the multiindex
% array (This may seem strange at first sight, but has many advantages in
% the implementation, very similar to that relation you have, when doing
% matrix multiplications).

fprintf([repmat('%7.3f ', 1, size(a_i_alpha,2)), '\n'], a_i_alpha)

%%
% The probability distribution of the germ is the joint distribution of the
% Semicircle and the Arcsine distributions. If we sample from the germ of
% V_a, we can see indeed that the marginal distributions are Semicircle and
% Arcsine, respectively.
xi = gpcgerm_sample(V_a, 30000);
plot_samples(xi); 
axis square;

%%
% We can also see that by looking at histograms of the marginals.
subplot(1,2,1); plot_density(xi(1,:), 'n', 30)
subplot(1,2,2); plot_density(xi(2,:), 'n', 30)

%%
% Now, we sample from the GPC expansion of the parameters, which gives the
% following joint distribution. To have a nicer plot with more evenly
% distributed points quasi Monte Carlo was specified by setting the 'mode'
% option to 'qmc'. (Note, that in the plot the arcsine distribution goes
% vertically with the larger spike on top, and the semicircle distribtion
% goes horizontally with the bump more to the left).
a_i_samples = gpc_sample(a_i_alpha, V_a, 30000, 'mode', 'qmc');
clf; plot_samples(a_i_samples);

%% Direct Monte Carlo sampling
% Now, as a first step, we can do a Monte-Carlo simulation of our
% deterministic model. Instead of making statistics, we'll just plot a
% bunch of samples, where the samples are this time created using Latin
% hypercube sampling (LHS). 

N = 30;
a_i_samples = gpc_sample(a_i_alpha, V_a, N, 'mode', 'lhs');
u_samples = zeros(N, size(pos,2));
minfo = diffusion_1d_init();
for i=1:N
    [u, minfo]=model_solve(minfo, a_i_samples(:,i));
    u_samples(i,:) = u;
end
plot(pos, u_samples)


%% Computing moments by sampling
% We have that also canned as a function. Showing mean and variance
% computed by MC and QMC.
[u_mean, u_var, minfo] = compute_moments_mc(minfo, a_i_alpha, V_a, 100);
subplot(1,2,1); plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
title('mc'); legend('mean-std', 'mean', 'mean+std'); ylim([0,3.5]); grid on;

[u_mean, u_var, minfo] = compute_moments_mc(minfo, a_i_alpha, V_a, 100, 'mode', 'qmc');
subplot(1,2,2); plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
title('qmc'); legend('mean-std', 'mean', 'mean+std'); ylim([0,3.5]); grid on;

%% Computing moments by quadrature
% Or we can compute that by projection/integration. With smolyak or tensor
% grid.
minfo = model_stats(minfo, 'reset');
[u_mean, u_var, minfo] = compute_moments_quad(minfo, a_i_alpha, V_a, 5, 'grid', 'smolyak');
model_stats(minfo, 'print');

subplot(1,2,1); 
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('smolyak'); ylim([0,3.5]); grid on;

minfo = model_stats(minfo, 'reset');
[u_mean, u_var, minfo] = compute_moments_quad(minfo, a_i_alpha, V_a, 5, 'grid', 'tensor');
model_stats(minfo, 'print');

subplot(1,2,2); 
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('tensor'); ylim([0,3.5]); grid on;

%% Response surfaces by projection
% First by projection
V_u = gpcbasis_create(V_a, 'p', 5);

minfo = model_stats(minfo, 'reset');
[u_i_alpha, minfo] = compute_response_surface_projection(minfo, a_i_alpha, V_a, V_u, 5, 'grid', 'full_tensor');
model_stats(minfo, 'print');

[u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
subplot(1,2,1); 
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('resp. surf. proj.'); ylim([0,3.5]); grid on;

subplot(1,2,2); 
plot_response_surface(u_i_alpha([10,30,60,90],:), V_u);
title('response surfaces at 0.1, 0.3, 0.6 and 0.9'); zlim([0, 5]);

%% Response surface by tensor grid interpolation
% Then by tensor grid interpolation
minfo = model_stats(minfo, 'reset');
[u_i_alpha, minfo] = compute_response_surface_tensor_interpolate(minfo, a_i_alpha, V_a, V_u, 5);
model_stats(minfo, 'print');

[u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
subplot(1,2,1); 
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('resp. surf. proj.'); ylim([0,3.5]); grid on;

subplot(1,2,2); 
plot_response_surface(u_i_alpha([10,30,60,90],:), V_u);
title('response surfaces at 0.1, 0.3, 0.6 and 0.9'); zlim([0, 5]);

%% Response surface by non-intrusive Galerkin
% Then by non-intrusive Galerkin (doesn't work)

max_iter = 20;
p_int = max(V_u{2}(:));
int_grid = 'full_tensor';
minfo.step_relax = 0.98;
minfo = model_stats(minfo, 'reset');
[u_i_alpha, minfo, x, w]=compute_response_surface_nonintrusive_galerkin(minfo, a_i_alpha, V_a, V_u, p_int, 'max_iter', max_iter, 'grid', int_grid);
model_stats(minfo, 'print_step_info');

[u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
subplot(1,2,1); 
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('resp. surf. proj.'); ylim([0,3.5]); grid on;

subplot(1,2,2); 
plot_response_surface(u_i_alpha([10,30,60,90],:), V_u);
title('response surfaces at 0.1, 0.3, 0.6 and 0.9'); zlim([0, 5]);

%% Response surface by intrusive Galerkin
% And now intrusive Stochastic Galerkin (just for the fun of it)
%
% First we need the stochastic Galerkin matrices for the two parameters
A_i = gpc_multiplication_matrices(a_i_alpha, V_a, V_u);
subplot(1,2,1); spy(A_i{1})
subplot(1,2,2); spy(A_i{2})

%%
% Then we can construct the full tensor product operator out of it.
K = tensor_operator_create({minfo.K, A_i});
clf; spy(tensor_operator_to_matrix(K))

%% 
% For applying stochastic Galerkin we need to apply the boundary conditions
% to the operator and right hand side. In sglib this is implemented such
% that the matrices (or tensor operators or whatever) are projected onto
% the inner (i.e. non-Dirichlet) nodes and the solution is later extended
% again by putting the values on the Dirichlet boundary back in.
P_I = minfo.P_I;
P_B = minfo.P_B;
N = size(minfo.K{1},1);
M = size(A_i{1},1);
Ki=apply_boundary_conditions_operator( K, P_I );


%%
% apply to the right hand side
G = zeros(N, M);
G(:,1) = minfo.g;
F = zeros(N, M);
F(:,1) = minfo.f;
Fi=apply_boundary_conditions_rhs( K, F, G, P_I, P_B );

%% 
% Here we transform the tensor operator to a matrix and solve with the
% standard matlab solve thing for sparse matrices
Ki_mat = tensor_operator_to_matrix(Ki);
Fi_vec = Fi(:);
Ui_vec = Ki_mat\Fi_vec;
Ui = reshape(Ui_vec, size(Fi));

%%
% Apply the boundary conditions and show the solution (and guess what? we
% get the same results as in the other cases)
u_i_alpha=apply_boundary_conditions_solution(Ui, G, P_I, P_B);

[u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
subplot(1,2,1); 
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('resp. surf. proj.'); ylim([0,3.5]); grid on;

subplot(1,2,2); 
plot_response_surface(u_i_alpha([10,30,60,90],:), V_u); 
title('response surfaces at 0.1, 0.3, 0.6 and 0.9'); zlim([0, 5]);


%%
% TODO: Compare response surface with true response
% TODO: Show grids for interpolation and projection
% TODO: Make plotting stuff into a function and explain
% TODO: Explain interpolation and/or projection shortly and show the code, then explain that that's been put into a function
% TODO: compare some samples computed directly and per surrogate model
% TODO: create model_stats(cmd) func (reset, print, ...)
% TODO: compare to dishis results
% TODO: create gpcbasis_info function (maybe remove gpcbasis_size)
