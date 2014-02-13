%%
clc
format compact; format long g
h=@normal_invcdf;
a=2;
f=@(x)(exp(abs(a*x)));

%exact result
2*exp(a^2/2)*normal_cdf(a)
% matlab's quad func
quad(@(x)(f(h(x))), 0, 1)
% gauss-hermite
gauss_hermite(f, 30)
% gauss-hermite again
[x,w]=polysys_int_rule('h', 30);
f(x)*w
% transformed gauss-legendre
[x,w]=polysys_int_rule('p', 300);
x=(x+1)/2; w=w;
f(h(x))*w
% trapezoidal geht gar nicht
[x,w]=trapezoidal_rule(20, 'interval', [0, 1]);
f(h(x))*w
% offene newton-cotes geht (resultat ist aber eher zufall)
[x,w]=newton_cotes_rule(16, 'open', true, 'interval', [0, 1]);
f(h(x))*w
% clenshaw-curtis does not work
[x,w]=clenshaw_curtis_rule(16, 'interval', [0, 1]);
f(h(x))*w
% fejer 1 should work
[x,w]=clenshaw_curtis_rule(30, 'mode', 'fejer1', 'interval', [0, 1]);
f(h(x))*w
% fejer 2 should work
[x,w]=clenshaw_curtis_rule(30, 'mode', 'fejer2', 'interval', [0, 1]);
f(h(x))*w


