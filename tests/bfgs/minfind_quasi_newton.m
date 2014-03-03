function [x,flag,iter]=minfind_quasi_newton(func, x0, varargin)

options=varargin2options(varargin,mfilename);
[maxiter, options]=get_option(options,'maxiter',100);
[abstol, options]=get_option(options,'abstol',1e-6);
[output_func, options]=get_option(options,'output_func',[]);
[verbosity, options]=get_option(options,'verbosity',0);
check_unsupported_options(options);


x=x0;
H = eye(length(x));

flag=1;
[y, dy] = funcall(func, x);
tol=abstol*max(1,tensor_norm(dy));
for iter=1:maxiter
    p = -H*dy;
    [alpha, xn, yn, dyn] = line_search_armijo(func, x, p, y, dy, 'alpha0', 0.6);
    %alpha = 1;
    %xn = x + alpha * p;
    if ~isempty(output_func)
        funcall(output_func, x, xn);
    end
    
    s = xn - x;
    yy = dyn - dy;
    [~, H] = qn_matrix_update('bfgs', [], H, yy, s);
    
    x = xn;
    y = yn;
    dy = dyn;
    if verbosity>0
        strvarexpand('minfind_quasi_newton: iter=$iter$, y=$y$, alpha=$alpha$ $eig(H)$');
    end
    if norm(dy)<tol
        flag=0;
        break;
    end
end    

