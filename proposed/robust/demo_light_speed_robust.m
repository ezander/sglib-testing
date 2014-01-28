function demo_light_speed_robust
% Generates some plots according to the Wikipedia article on robust
% statistics http://en.wikipedia.org/wiki/Robust_statistics

global seed
seed = rand('seed');

% Estimation of Location
% http://en.wikipedia.org/wiki/Robust_statistics#Estimation_of_location

time = read_data('light');
[x,p]=kernel_density(time);

subplot(2,2,1)
plot(x,p);
hold all;
rug_plot(time);
hold off;


subplot(2,2,2)
y = sort(time);
n = length(time);
x = normal_invcdf(linspace(0.5/n, 1 - 0.5/n, n));
plot(x, y, 'x');

subplot(2,2,3)
mu = bootstrap_mean(time, 10000);
kernel_density(mu);

subplot(2,2,4)
mu = bootstrap_tmean(time, 10000);
kernel_density(mu);



function mu=bootstrap_mean(x, m)
rnginit();
y=bootstrap_sample(x, m);
mu=mean(y,2);

function mu=bootstrap_tmean(x, m)
rnginit();
y=bootstrap_sample(x, m);
y=sort(y,2);
n = length(x);
n10 = round(n*0.1);
mu=mean(y(:,1+n10:n-n10),2);

function rnginit()
global seed
disp(seed)
rand('seed', seed);
