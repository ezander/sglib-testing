function test_optim

func = @rosenbrock;
x0=[-1.9;2];
%x0=[4*rand-2;4*rand-2];
x_ex = [1; 1];

clc
plot_func(func);
outfunc = funcreate(@iterplot, @funarg, @funarg, 'ro-');
[x,flag,iter] = minfind_newton(func, x0, 'abstol', 1e-8, 'output_func', outfunc, 'verbosity', 1)
outfunc = funcreate(@iterplot, @funarg, @funarg, 'go-');
[x,flag,iter] = minfind_quasi_newton(func, x0, 'abstol', 1e-8, 'output_func', outfunc, 'verbosity', 1)



function plot_func(func)
N=100;
x=linspace(-2,2,N);
y=linspace(-1,3,N);
x=linspace(-2,2,N);
y=linspace(-8,5,N);
[X,Y]=meshgrid(x,y);
z = funcall(func, [X(:),Y(:)]');
Z=reshape(z, size(X));
hold off
%plotedit on
contour(X,Y,Z,[1:10].^2);
