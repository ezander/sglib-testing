function test_tensor_optim

clc
format compact; format short g;

rand_seed(12345)

N=3;
M=4;
K=5;
A=TensorOperator({spdmat(N), spdmat(M), spdmat(K); spdmat(N), spdmat(M), spdmat(K)});
x_ex=rand(N,M,K);
b=A*x_ex;

% solve by converting to full system
x=reshape(A.asmatrix()\tensor_to_vector(b), tensor_size(x_ex));
tensor_error(x, x_ex, 'relerr', true)

eig(A.asmatrix())
cond(A.asmatrix())

% solve by quasi newton stuff
qnewton_opts.abstol = 1e-5;
qnewton_opts.verbosity = 1;
qnewton_opts.maxiter= 1000;
if false
    func = funcreate(@myfunctional, A, b, @funarg);
    qnewton_opts.line_search_func = @line_search_armijo;
    qnewton_opts.line_search_opts = {'stretch', true, 'alpha0', 1.5};
else
    func = funcreate(@myfunctional_grad, A, b, @funarg);
    qnewton_opts.line_search_func = @line_search_gradient;
    qnewton_opts.line_search_opts = {'alpha0', 2};
end
[func, stats] = CallStatistics.wrap_function(func);
I=IdentityOperator.from_vector(b);
H0 = inv(LBFGSOperator(I, 30));
x0 = tensor_null(b);
x = minfind_quasi_newton(func, x0, H0, qnewton_opts);
stats
tensor_error(x, x_ex, 'relerr', true)


function [f,df]=myfunctional_grad(A, b, x)
df = -tensor_add(b, A*x, -1);
f = nan;

function [f,df]=myfunctional(A, b, x)
y=A*x;
r = -tensor_add(b, y, -1);
f = 0.5*tensor_scalar_product(x,y) - tensor_scalar_product(x,b);
df=r;


function A=spdmat(n)
A=rand(n);
A=A'*A+0.1*eye(n);
