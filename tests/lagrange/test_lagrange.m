n=8;

fn = rand(n,1);
xn = [0; rand(n-2,1); 1];
%newton_tableau(x, f)
x = linspace(0,1);

y=lagrange_newton(xn, fn, x)
plot(x,y, xn, fn, 'rx');

