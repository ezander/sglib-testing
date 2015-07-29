function test_bayes_mcmc
% TEST_BAYES_MCMC Test code for the Bayes MCMC sampling method

prior_dist=gendist_create('normal', {0, 2});
error_dist=gendist_create('normal', {0, 0.1});
[~,sig2]=gendist_moments(prior_dist);
prop_dist = gendist_create('uniform', {-sqrt(sig2), sqrt(sig2)});

N=10000;
y_m=1;
x0 = gendist_moments(prior_dist);
%func=@(x)(gendist_pdf(x, prior_dist));
func=@(x)(gendist_pdf(y_m-x.^2, error_dist).*gendist_pdf(x, prior_dist));

%X=mh_sample(N, x0, func, prop_dist);
X=mh_sample_parallel(N, x0, func, prop_dist);
X_prior=gendist_sample(N, prior_dist);


multiplot_init(1,1)

% Plot MH density and correlation
multiplot;
hold off;
plot_density(X); hold all; grid on;
legend_add('metropolis');
plot_density(prior_dist);
legend_add('exact dist');


function X=mh_sample(N, x0, func, prop_dist)
% MH_SAMPLE Basic version of the Metropolis-Hastings sampler
M = 1000;
x=x0;
X=[];
for i=1:N+M
    xn=x+gendist_sample(1,prop_dist);
    a=funcall(func,xn)/funcall(func,x);
    a=a*gendist_pdf(x-xn, prop_dist)/gendist_pdf(xn-x,prop_dist);
    if a>=1 || rand<a
        x=xn;
    end
    if i>M
        X=[X x];
    end
end
X=X(1:N);


function X=mh_sample_parallel(N, x0, func, prop_dist)
% MH_SAMPLE_PARALLEL Parallel version of the Metropolis-Hastings sampler

K = 100;
M = K*ceil(1000/K);

x = repmat(x0, 1, K);
X=[];
p = funcall(func, x);

for i=1:ceil((N+M)/K)
    xn=x+gendist_sample(K, prop_dist)';
    pn=funcall(func, xn);
    a=pn./p;
    a=a.*gendist_pdf(x-xn, prop_dist)./gendist_pdf(xn-x,prop_dist);
    ind = (a>=1 | rand(1,K)<a);
    x(ind)=xn(ind);
    p(ind)=pn(ind);
    if i>M/K
        X=[X, x];
    end
end
X=X(1:N);

