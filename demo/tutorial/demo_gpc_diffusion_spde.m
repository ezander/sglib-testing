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
plot(pos, u); enhance_plot;

%% Stochastic description of the parameters
% Now suppose the parameters a1 and a2 are uncertain and described by
% by random variables. To make things a bit interesting, we take two
% different distributions for a1 and a2, with very different shapes and
% very different ranges.
% 
% Suppose therefore, that a1 follows a Beta(1.2,2) distribution, which is
% shifted and scaled such that it has support [0.5, 5]. In sglib this can be
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

plot_density(a1_dist, 'type', 'both');
legend( 'pdf a1', 'cdf a1'); enhance_plot;

%%
% The second parameter a2 follows a shifted and scaled Beta(0.7, 0.6)
% distribution with support on [50, 150]. The distribution of a2 is also
% asymetrical, but it is bimodal and has sharp singularities at both ends.
% In order to better see the singularities in the plot some extra x-values
% to plot in their vicinity are specified with the 'extra_x' option.

a2_dist = gendist_create('beta', {0.7, 0.6});
a2_dist = gendist_fix_bounds(a2_dist, 50, 150);

plot_density(a2_dist, 'type', 'both', 'extra_x', [50.0001, 149.999])
legend( 'pdf a2', 'cdf a2'); enhance_plot;

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
enhance_plot;

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
% The moment approximations looked quite good. Now, we plot a kernel
% density estimate of the GPC approximation of a2 and compare that to the
% true distribution. Note that for the kernel density estimation we need to
% specify the with of the kernels (setting 'kde_sig' to 0.5), because the
% default used here (Silverman's rule of thumb) is not appropriate for this
% kind of bimodal distribution.
a2_samples = gpc_sample(a2_alpha, V2, 100000);

plot_density(a2_samples, 'type', 'kernel', 'kde_sig', 0.5, 'rug', true, 'max_rug', 300);
plot_density(a2_dist, 'hold', true);
legend('GPC kernel density', 'GPC samples', 'Exact pdf'); enhance_plot;


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
plot_density(xi(1,:), 'n', 30); enhance_plot;
plot_density(xi(2,:), 'n', 30); enhance_plot;

%%
% Now, we sample from the GPC expansion of the parameters, which gives the
% following joint distribution. To have a nicer plot with more evenly
% distributed points quasi Monte Carlo was specified by setting the 'mode'
% option to 'qmc'. (Note, that in the plot the arcsine distribution goes
% vertically with the larger spike on top, and the semicircle distribtion
% goes horizontally with the bump more to the left).
a_i_samples = gpc_sample(a_i_alpha, V_a, 30000, 'mode', 'qmc');
plot_samples(a_i_samples); enhance_plot;

%% Direct Monte Carlo sampling
% Now, as a first step, we can do a Monte-Carlo simulation of our
% deterministic model. Instead of making statistics, we'll just plot a
% bunch of samples, where the samples are this time created using Latin
% hypercube sampling (LHS). 

N = 30;
a_i_samples = gpc_sample(a_i_alpha, V_a, N, 'mode', 'lhs');
u_samples = zeros(N, size(pos,2));
model = diffusion_1d_init();
for i=1:N
    [u, model]=model_solve(model, a_i_samples(:,i));
    u_samples(i,:) = u;
end
plot(pos, u_samples); enhance_plot;

%%
% Note: the function diffusion_1d_init sets up a model structure, which contains
% information about the deterministic model, e.g. its dimension, how many parameters it
% has, a handle to the solve function (for MC and collocation), a handle to
% a solver step function (for non-intrusive Galerkin), etc. The advantage
% is that the complete information about the model is captured in one
% structure and only this struct needs to be passed around.

%% Computing moments by sampling
% Monte Carlo is can also be invoked directly for the model in form of some
% ready made function. The following shows how to compute the mean and
% variance by MC and QMC.
[u_mean, u_var, model] = compute_moments_mc(model, a_i_alpha, V_a, 100);
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
title('mc'); legend('mean-std', 'mean', 'mean+std'); ylim([0,3.5]); grid on;
enhance_plot;

[u_mean, u_var, model] = compute_moments_mc(model, a_i_alpha, V_a, 100, 'mode', 'qmc');
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
title('qmc'); legend('mean-std', 'mean', 'mean+std'); ylim([0,3.5]); grid on;
enhance_plot;

%% Computing moments by quadrature
% There are also function that do that by by quadrature. The following code
% shows how to do that with a Smolyak and with a tensor grid.
model = model_stats(model, 'reset');
[u_mean, u_var, model] = compute_moments_quad(model, a_i_alpha, V_a, 5, 'grid', 'smolyak');
model_stats(model, 'print');

plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('smolyak'); ylim([0,3.5]); grid on; enhance_plot;

model = model_stats(model, 'reset');
[u_mean, u_var, model] = compute_moments_quad(model, a_i_alpha, V_a, 5, 'grid', 'tensor');
model_stats(model, 'print');

plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('tensor'); ylim([0,3.5]); grid on; enhance_plot;


%% Response surfaces by projection
% While the previous functions directly computed moments of the models
% response it is also possible to compute response surfaces. The code first
% constructs a function space V_u for the response surface (or surrogate
% model if you prefer) of orthogonal polynomials defined on the same germ
% as a but now up to complete degree 5. Then the response surface is
% computed.
V_u = gpcbasis_create(V_a, 'p', 5);

model = model_stats(model, 'reset');
[u_i_alpha, model] = compute_response_surface_projection(model, a_i_alpha, V_a, V_u, 5, 'grid', 'smolyak');
model_stats(model, 'print');

%%
% To visualise the result first the mean and variance are plotted. And then
% the response surfaces at the points x=0.1, 0.3, 0.6, and 0.9 are plotted.
% Since the parameter space is two dimensional those are indeed surfaces.
[u_mean, u_var] = gpc_moments(u_i_alpha, V_u);
plot(pos, u_mean-sqrt(u_var), pos, u_mean, pos, u_mean+sqrt(u_var));
legend('mean-std', 'mean', 'mean+std'); 
title('resp. surf. proj.'); ylim([0,3.5]); grid on; enhance_plot;

plot_response_surface(u_i_alpha([10,30,60,90],:), V_u);
title('response surfaces at 0.1, 0.3, 0.6 and 0.9'); zlim([0, 5]); enhance_plot;

%% Response surface by tensor grid interpolation
% Here the response surfaces are computed by tensor grid interpolation. In
% order to reduce the code the plotting commands have been moved into a
% script, since they will be the same for all response surface
% computations.
model = model_stats(model, 'reset');
[u_i_alpha, model] = compute_response_surface_tensor_interpolate(model, a_i_alpha, V_a, V_u, 5);
model_stats(model, 'print');

plot_response_surface_results

%% Response surface by non-intrusive Galerkin
% Here the response surfaces are computed via the non-intrusive Galerkin
% method.
p_int = max(V_u{2}(:))+1;
int_grid = 'full_tensor';
[u_i_alpha, model, x, w]=compute_response_surface_nonintrusive_galerkin(model, a_i_alpha, V_a, V_u, p_int, 'maxiter', maxiter, 'grid', int_grid);

plot_response_surface_results


%% Response surface by intrusive Galerkin
% The computation of the response surfaces by intrusive stochastic
% Galerking is a bit more complicated and will be shown here without the
% accompanying math formulas. 
%
% First we need the stochastic Galerkin matrices for the two parameters a1
% and a2, 
A_i = gpc_multiplication_matrices(a_i_alpha, V_a, V_u);
subplot(1,2,1); spy(A_i{1})
subplot(1,2,2); spy(A_i{2})

%%
% Then we can construct the full tensor product operator out of it.
K = tensor_operator_create({model.K, A_i});
clf; spy(tensor_operator_to_matrix(K))

%%
% Setup the boundary conditions and right hand side (very manually here,
% since they are deterministic, for stochastic BCs and RHSs are some
% specialised functions available).
[N,M]=operator_size(K, 'domain_only', true, 'contract', false);
G = zeros(N, M);
G(:,1) = model.g;
F = zeros(N, M);
F(:,1) = model.f;

%% 
% For applying stochastic Galerkin we need to apply the boundary conditions
% to the operator and right hand side to get a modified system with
% boundary conditions included. 

P_I = model.P_I;
P_B = model.P_B;
[Kn,Fn]=apply_boundary_conditions_system( K, F, G, P_I, P_B );

%% 
% There are different ways now to solve this linear system. First note that
% the operator is still in tensor product form

Kn

% Here we transform the tensor operator to a matrix and solve with the
% standard matlab solve thing for sparse matrices
Kn_mat = tensor_operator_to_matrix(Kn);
Fn_vec = Fn(:);
U_vec = Kn_mat\Fn_vec;
u_i_alpha=reshape(U_vec, size(F));
plot_response_surface_results

%%
% Here we construct a preconditioner for later use in a PCG solver. (This
% would be easier if we started with a KL expansion of the coefficient
% field a, but now it requires bit handiwork).
a1_mean = gendist_moments(a1_dist);
a2_mean = gendist_moments(a2_dist);
P = a1_mean * K{1,1} + a2_mean*K{2,1};
Pn=apply_boundary_conditions_system(P, F, G, P_I, P_B);
Pinv = operator_from_matrix_solve(Pn);

%%
% Now solve with a PCG and after that with a solver based on simple
% iterations (or preconditioned Richardson if you prefer). The response
% surfaces come out identical.
[Un, flag, info] = tensor_solve_pcg(Kn, Fn, 'Minv', Pinv);
tensor_solver_message(info);

u_i_alpha=Un;
plot_response_surface_results


[Un, flag, info] = tensor_solve_simple(Kn, Fn, 'Minv', Pinv);
tensor_solver_message(info);

%%
% Here an standard Matlab iterative solver is used (@pcg) for solving the
% system. Since Matlab's solver can't handle the sglib data structures and
% operators directly tensor_solve_matlab_wrapper constructs a wrappers and
% passed them to the Matlab pcg.
[Pinv, P] = stochastic_preconditioner(Kn, 'precond_type', 'vanloan', 'num_iter', 5);
Un = tensor_solve_matlab_wrapper(@pcg, Kn, Fn, 'Minv', Pinv);
u_i_alpha=Un;
plot_response_surface_results
