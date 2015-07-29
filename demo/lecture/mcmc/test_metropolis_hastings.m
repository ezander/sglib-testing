function test_metropolis_hastings
% TEST_METROPOLIS_HASTINGS Test code for the Metropolis-Hastings sampling method
% See: https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm

%dist=gendist_create('normal', {20, 3});
dist=gendist_create('beta', {2, 4});
%dist=gendist_create('lognormal', {-1, .37});

[~,sig2]=gendist_moments(dist);
prop_dist = gendist_create('normal', {-0.02, 0.3*sqrt(sig2)});

N=10000;

%X=mh_sample(N, dist, prop_dist);
X=mh_sample_parallel(N, dist, prop_dist);
X2=gendist_sample(N, dist);


multiplot_init(2,2)

% Plot MH density and correlation
multiplot;
hold off;
plot_density(X); hold all; grid on;
legend_add('metropolis');
plot_density(dist);
legend_add('exact dist');

multiplot;
plot(X(1:end-1), X(2:end), 'x');


% Plot density and correlation from direct sampling
multiplot;
plot_density(X2); hold all; grid on;
legend_add('direct samples');
plot_density(dist);
legend_add('exact dist');

multiplot;
plot(X2(1:end-1), X2(2:end), 'x');

multiplot_adjust_range('separate', 'rows');



function X=mh_sample(N, dist, prop_dist)
% MH_SAMPLE Basic version of the Metropolis-Hastings sampler
M = 1000;
x=gendist_moments(dist);
X=[];
for i=1:N+M
    xn=x+gendist_sample(1,prop_dist);
    a=gendist_pdf(xn,dist)/gendist_pdf(x,dist);
    a=a*gendist_pdf(x-xn, prop_dist)/gendist_pdf(xn-x,prop_dist);
    if a>=1 || rand<a
        x=xn;
    end
    if i>M
        X=[X x];
    end
end
X=X(1:N);


function X=mh_sample_parallel(N, dist, prop_dist)
% MH_SAMPLE_PARALLEL Parallel version of the Metropolis-Hastings sampler

K = 100;
M = K*ceil(1000/K);

x=gendist_moments(dist);
x = repmat(x, 1, K);
X=[];
p = gendist_pdf(x, dist);

for i=1:ceil((N+M)/K)
    xn=x+gendist_sample(K, prop_dist)';
    pn=gendist_pdf(xn, dist);
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

