function [L, n]=covariance_decomp(C, varargin)

options = varargin2options(varargin, mfilename);
[fillup, options]=get_option(options, 'fillup', false);
[always_lower, options]=get_option(options, 'always_lower', false);
check_unsupported_options(options);

[L,p]=chol(C, 'lower');
if p>0
    % that indicates chol did not work 
    
    if true
        [U,D]=eig(C);
    else
        [U,D]=ldl(C);
    end
    d = diag(D);
    min_d = abs(max(d)) * 1e-14;
    ind = (d>min_d);
    L = binfun(@times, U(:,ind), sqrt(d(ind)'));
    
    if always_lower
        % Note: is X=qr(A), then triu(X) is the upper triangular part of 
        % R where Q*R=A is the QR decomposition of A
        L=triu(qr(L'))';
    end
end

[m, n] = size(L);
if fillup && n<m
    L = [L, zeros(m, m-n)];
end

assert( norm(L*L'-C)<1e-12*norm(C) )
