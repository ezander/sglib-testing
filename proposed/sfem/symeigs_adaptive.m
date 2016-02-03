function [V, d] = symeigs_adaptive(A, B, varargin)



n = size(A,1);

options = varargin2options(varargin, mfilename);
[eig_ratio, options]=get_option(options, 'eig_ratio', 0);
[eigs_disp, options]=get_option(options, 'eigs_disp', 0);
[sum_eigs, options]=get_option(options, 'sum_eigs', inf);
[sum_sq_eigs, options]=get_option(options, 'sum_sq_eigs', inf);
[num_eigs, options]=get_option(options, 'num_eigs', n);
[stop_fun, options]=get_option(options, 'stop_fun', []);
[k_min, options]=get_option(options, 'k_min', 40);
[k_max, options]=get_option(options, 'k_max', 400);
[k_start, options]=get_option(options, 'k_start', k_min);
[num_pred, options]=get_option(options, 'num_pred', k_min/2);
check_unsupported_options(options);

opts = struct();

V = zeros(n,0);
d = zeros(0,1);
k1 = inf;
Adef = A;
k = k_min;
stats.k = [];
while true
    opts.disp = eigs_disp;
    stats.k(end+1) = k;
    stats.k
    if ~false
        [Vn, Dn] = eigs(Adef, B, k, 'lm', opts);
        dn = diag(Dn);
        d = [d; dn]; %#ok<AGROW>
        V = [V, Vn]; %#ok<AGROW>
        M = diag(sqrt(dn)) * Vn' * B;
        Adef = Adef - M' * M;
        assert(all(all(Adef==Adef')));
    else
        [V, D] = eigs(Adef, B, sum(stats.k), 'lm', opts);
        d = diag(D);
    end
        
    
    
    
    %%
    dl = d(end-num_pred+1:end);
    [a, b] = polyfit(-num_pred:-1, log(dl'), 1);
    alpha = 0.5;
    dp = [d; exp(alpha*a(1)*(0:k_max)'+a(2))];
    semilogy(dp);
    
    %%
    
    k1 = inf;
    if isfinite(sum_eigs) && sum(abs(dp))>sum_eigs
        k0 = find(cumsum(abs(dp))>=sum_eigs, 1, 'first');
        k1 = min(k1, k0);
    end
    
    if isfinite(sum_sq_eigs) && sum(dp.^2)>sum_sq_eigs
        k0 = find(cumsum(dp.^2)>=sum_sq_eigs, 1, 'first');
        k1 = min(k1, k0);
    end

    if dp(end) <= eig_ratio*d(1)
        k0 = find(dp<=eig_ratio*d(1), 1, 'first');
        k1 = min(k1, k0);
    end
    
    if length(dp) >= num_eigs
        k0 = num_eigs;
        k1 = min(k1, k0);
    end
    
    if ~isempty(stop_fun)
        k0 = funcall(stop_fun, dp);
        k1 = min(k1, k0);
    end
    
    if isfinite(k1)
        if k1<=length(d)
            break;
        end
        k = min(num_eigs, k1) - length(d);
        k = min(max(k, k_min), k_max);
    else
        k = k_max;
    end
end
V = V(:, 1:k1);
d = d(1:k1);
strvarexpand( 'k1: $k1$, sum: $sum(stats.k)$, array: $stats.k$');
