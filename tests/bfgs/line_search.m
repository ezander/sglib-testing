function [alpha, xn, yn, dyn, flag]=line_search(func, x, p, y, dy, varargin)
% LINE_SEARCH Performs a ... line search.
%   [ALPHA,XN,YN,DYN,FLAG]=LINE_SEARCH(FUNC,X,P,Y,DY,OPTIONS)
%
%   The algorithm works by ... (see [1], Algorithm 3.2, page 59) for the
%   line search algorithm and (see [1], Algorithm 3.3, page 60).
%
% Options:
%   alpha0: double {1}
%     Initial value for alpha.
%   rho: double {0.5}
%     Reduction factor for alpha.
%   c: double {1e-4}
%     Constant in the Armijo condition.
%   maxiter: integer {100}
%     Maximum number of iterations until suitable alpha is found.
%   verbosity: integer {0}
%     If larger zero diagnoistic messages will be printed.
%   
% References
%   [1] Nocedal, J. and Wright, S.J. (1999): Numerical Optimization,
%       Springer-Verlag. ISBN 0-387-98793-2.


options=varargin2options(varargin,mfilename);
[alpha0, options]=get_option(options,'alpha0',1);
[alpha_max, options]=get_option(options,'alpha_max',10);
[rho, options]=get_option(options,'rho',0.5);
[c1, options]=get_option(options,'c1',1e-4);
[c2, options]=get_option(options,'c2',0.9);
[maxiter, options]=get_option(options,'maxiter',100);
[verbosity, options]=get_option(options,'verbosity',0);
check_unsupported_options(options);

% parameters
alpha1 = (pi)/8;

% optional stuff
alpha0 = 0;

clc
format compact
format short g
[y0, dy0] = funcall(phi_func, alpha0)

for i=1:maxiter
    [yi, dyi] = funcall(phi_func, alpha1)
    alpha1
    
    % No sufficient decrease (would be yi<=y0 + x1 * alpha1 * dy0, i.e. the
    % Armijo condition) or an increase with respect to the last iterate
    % Then we nees to zoom in as in between 0 and alpha1 there must be a
    % point where the conditions are fulfilled.
    if yi>(y0 + c1 * alpha1 * dy0) || (i>1 && yi>=yo)
        alpha = zoom(alpha0, alpha1);
        break;
    end
    
    % Sufficient decrease is there, now check the curvature condition. If
    % that's fulfilled we can exit.
    if abs(dyi) <= c2 * abs(dy0)
        alpha = alpha1;
        break;
    end
    
    % If we have a positive derivative we need to zoom in, too.
    if dyi>=0 
        alpha = zoom(alpha1, alpha0);
        break;
    end
    
    % Better choice for alpha1?
    alpha0 = alpha1;
    alpha1 = min(2 * alpha1, 0.5 * (alpha1 + alpha_max));
    yo = yi;
end

function alpha = zoom(alpha1, alpha2)
alpha = 0.5 * (alpha1 + alpha2);


function [y,dy] = test_func(x)
y = -sin(x);
dy = -cos(x);
