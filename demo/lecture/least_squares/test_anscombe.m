function test_anscombe

%show_stats
%show_fits(1)
%show_errors(3, 12)
%show_l1_fits(false)

function show_errors(num, max_n)
for n=1:max_n
    [x,y] = anscombes_quartet(num);
    q = polynomial_regression(x, y, n);
    [n, evaluate_mse(x, y, q)]
end

function show_fits(n)
multiplot_init(4);
for num = 1:4
    multiplot
    [x,y] = anscombes_quartet(num);

    scatter(x,y)
    q = polynomial_regression(x, y, n);
    plot_poly_curve(q, min(x), max(x));
    %p
    evaluate_mse(x, y, q)
end
%multiplot_adjust_range()

function show_l1_fits(init)
if init
    multiplot_init(4);
end
for num = 1:4
    multiplot
    [x,y] = anscombes_quartet(num);
    if ~init; hold all; end
    scatter(x,y)


    %%
    [Q0, Q1] = meshgrid(linspace(0,10), linspace(0,1));
    res = reshape(sum(abs(Q0(:) + Q1(:)*x' - y'),2), size(Q0));
    %surf(Q0, Q1, res)
    %shading interp
    %grid on
    %%
    [a,b]=min(res(:));
    q = [Q1(b); Q0(b)];
    plot_poly_curve(q, min(x), max(x), 'color', 'g');
    %p
    evaluate_mse(x, y, q)
end
%multiplot_adjust_range()



function show_stats
n=1;
for num = 1:4
    [x,y] = anscombes_quartet(num);
    q = polynomial_regression(x, y, n);
    
    mu_y = mean(y);
    ss_tot = sum((y - mu_y).^2);
    ss_exp = sum((mu_y - polyval(q,x)).^2);
    ss_res = sum((y - polyval(q,x)).^2);
    R2 = 1 - ss_res / ss_tot;
    %[ss_exp, ss_res, ss_exp+ss_res, ss_tot]
    [mean(x), var(x), round(mean(y), 2), round(var(y), 3), round(corr(x,y),3), round(q(1:2)',2), round(R2,2)]
end



function q = polynomial_regression(x, y, n)
A = poly_matrix(x, n);
s=warning('off', 'MATLAB:rankDeficientMatrix');
q = A\y;
warning(s);

function mse=evaluate_mse(x, y, q)
yq = polyval(q, x);
mse=sum( (y-yq).^2 );
