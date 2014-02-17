function test_ziggurat
dist=gendist_create('normal', {0,1});
xi1=gendist_sample(1000, dist);
xi2=sample_zigg(dist, 1000000);
plot_density(xi2, 30)


function xi=sample_zigg(dist, N)

%%
%[xb, xn, yb, yn] = get_ziggurat_slices(dist, 3.4)
[xb, xn, yb, yn] = get_ziggurat_slices(dist, 1.65);
xi = sample_from_ziggurat(dist, N, xb, xn, yb, yn);

%%
function xi=sample_from_ziggurat(dist, N, xb, xn, yb, yn)
nl = length(xn);
xi = nan(1,N);
Nr = sum(isnan(xi));
while Nr>0
    n = randi(nl, 1, Nr);
    u0 = rand(1, Nr);
    xi0 = u0 .* xb(n);
    
    ind_ok = (xi0<=xn(n));
    xi0(~ind_ok)=nan;
    xi(isnan(xi)) = xi0;
    n = n(~ind_ok);
    u0 = u0(~ind_ok);
    xi0 = u0 .* xb(n);
    Nr = sum(isnan(xi));
    if ~Nr; break; end
    
    u1 = rand(1, Nr);
    y =  yb(n) + u1 .* (yn(n) - yb(n));
    ind_ok = (y < gendist_pdf(xi0, dist));

    xi0(~ind_ok)=nan;
    xi(isnan(xi)) = xi0;
    n = n(~ind_ok);
    u0 = u0(~ind_ok);
    y = y(~ind_ok);
    xi0 = u0 .* xb(n);
    Nr = sum(isnan(xi));

    xi0(:) = nan;
    xi0(n==1) = tail_sample(dist, xb(1), sum(n==1));
    xi(isnan(xi)) = xi0;
    
    % need to make the tail thing here
    Nr = sum(isnan(xi));
end    

%%


function [xb, xn, yb, yn] = get_ziggurat_slices(dist, x)
y = gendist_pdf(x, dist);
A = x * gendist_pdf(x, dist) + (gendist_cdf(inf, dist) - gendist_cdf(x, dist));
xb = [A/y];
xn = [x];
yb = [0];
yn = [y];
while true
    xb(end+1)=x;
    yb(end+1) = y;
    dy = A / x;
    y = y + dy;
    yn(end+1) = y;
    if y > gendist_pdf(0,dist)
        xn(end+1) = 0;
        break;
    end
    x = fzero(@(x)(gendist_pdf(x,dist)-y), [0, x]);
    xn(end+1) = x;
end


function xi = tail_sample(dist, from_x, n)
min_y = gendist_cdf(from_x, dist);
u = min_y + (1-min_y)*rand(n,1);
xi=gendist_invcdf(u, dist);


