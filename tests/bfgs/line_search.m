function alpha = line_search()

% p 59 algorithm 3.2

% parameters
alpha1 = (pi)/8;
alpha_max = 10;
phi_func = @test_func; % maybe from vec + func

% optional stuff
c1 = 1e-4;
c2 = 0.9;
maxiter=100;


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
