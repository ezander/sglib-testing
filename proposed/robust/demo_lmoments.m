% Computes some L-moments and L-ratios according to the Wikipedia article
% http://en.wikipedia.org/wiki/L-moment 

% Set some parameters (formatting and number of random samples used)
format_str = '%s: l1 =% 6.4f, l2 =% 6.4f, tau3 =% 6.4f, tau4 =% 6.4f\n';
N = 10000;


% Uniform
a = 0; b = 1;
x = (b-a) * rand(1, N) + a;
lm_uniform_exact = [(a+b)/2, (b-a)/6, 0, 0];
lm_uniform_empiric = lmoments(x);

disp(' ');
fprintf(format_str, 'Uniform (exact)  ', lm_uniform_exact)
fprintf(format_str, 'Uniform (empiric)', lm_uniform_empiric)


% Normal
mu = 0; sig2 = 1;
x = randn(1, N) * sqrt(sig2) + mu;
lm_normal_exact = [mu, sqrt(sig2/pi), 0, 0.1226];
lm_normal_empiric = lmoments(x);

disp(' ');
fprintf(format_str, 'Normal (exact)  ', lm_normal_exact)
fprintf(format_str, 'Normal (empiric)', lm_normal_empiric)


% Exponential
lambda=1.3;
x = exponential_invcdf(rand(1, N), lambda);
lm_exponential_exact = [1/lambda, 0.5/lambda, 1/3, 1/6];
lm_exponential_empiric = lmoments(x);

disp(' ');
fprintf(format_str, 'Exponential (exact)  ', lm_exponential_exact)
fprintf(format_str, 'Exponential (empiric)', lm_exponential_empiric)
