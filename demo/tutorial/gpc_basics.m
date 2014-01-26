%% Creating GPC bases and integration

%% Create a rather large basis
% The following creates a GPC basis in five variable from which the first
% is Hermite/Gaussian, the second is Uniform/Legendre, the third is again
% Hermite/Gaussian but with normalised Hermite polynomials, the fourth is
% Laguerre/Exponential and the fifth is ChebyshevT/Arcsine. 
V = gpcbasis_create('HphLT', 'p', 40);
%%
% The letters follow the standard symbol of the family of orthogonal polynomials, where a
% capital letter denotes the standard definition, and a small letter the
% normalised version of the polynomials. The "'p', 40" tells the functions
% to use total polynomials up to degree 40.

%% 
% Lets check how many basis functions we have:
fprintf('%d basis functions in %d RVs\n', ...
    gpcbasis_size(V,1), gpcbasis_size(V,2));

%%
% We can integrate a function over this space. The integration is done with
% a tensorised Gauss rule (i.e. Gauss-Hermite x Gauss-Legendre x
% Gauss-Hermite x Gauss-Laguerre x Gauss-Chebyshev) on a Smolyak grid and
% on a full tensor grid.
func = @(x)(100*cos(x(1,:)).*cos(x(2,:)).*cos(x(3,:)).*sin(x(4,:)).*cos(x(5,:)));
for p=1:8
    Is=gpc_integrate(func, V, p);
    It=gpc_integrate(func, V, p, 'grid', 'tensor');
    fprintf('p=%2d I(smolyak)=% +6.2f I(tensor)=% +6.2f\n', p, Is, It);
end

%%
% We can do the same with Monte Carlo
N=100000;
xi = gpcgerm_sample(V, N);
Imc = mean(funcall(func, xi));
fprintf('I_mc=% +6.2f \n', Imc);

%%
% or Latin hypercube
xi = gpcgerm_sample(V, N, 'mode', 'lhs');
Ilhs = mean(funcall(func, xi));
fprintf('I_lhs=% +6.2f \n', Ilhs);

%%
% or with quasi Monte Carlo (Halton)
xi = gpcgerm_sample(V, N, 'mode', 'qmc');
Iqmc = mean(funcall(func, xi));
fprintf('I_qmc=% +6.2f \n', Iqmc);


%% A smaller basis for which its easier to plot something
% First show the multiindices of gpc basis function 
V=gpcbasis_create('Hp', 'p', 4);
I=V{2}

%%
% Sampling from this basis (you see that its Gaussian in x-direction
% uniform in y-direction)
xi = gpcgerm_sample(V, 100000);
plot(xi(1,:), xi(2,:), '.');

%%
% Show a smolyak grid for the combination
[x,w]=gpc_integrate([], V, 5);
plot(x(1,:), x(2,:), 'x');

%%
% And now a full tensor grid for ChyshevU/Laguerre with orders 3 and 8
% respectively
[x,w]=gpc_integrate([], gpcbasis_create('ul'), [3, 8], 'grid', 'tensor');
plot(x(1,:), x(2,:), 'x');

%% foo
