function test_anscombe

show_stats
%show_fits(8)
%show_errors(4, 12)


function show_errors(num, max_n)
for n=1:max_n
    [x,y] = anscombes_quartet(num);
    p = polynomial_regression(x, y, n);
    [n, evaluate_mse(x, y, p)]
end

function show_fits(n)
multiplot_init(4);
for num = 1:4
    multiplot
    [x,y] = anscombes_quartet(num);

    scatter(x,y)
    p = polynomial_regression(x, y, n);
    plot_poly_curve(p, min(x), max(x));
    %p
    evaluate_mse(x, y, p)
end
%multiplot_adjust_range()

function show_stats
n=1;
for num = 1:4
    [x,y] = anscombes_quartet(num);
    p = polynomial_regression(x, y, n);
    
    mu_y = mean(y);
    ss_tot = sum((y - mu_y).^2);
    ss_exp = sum((mu_y - polyval(p,x)).^2);
    ss_res = sum((y - polyval(p,x)).^2);
    R2 = 1 - ss_res / ss_tot;
    %[ss_exp, ss_res, ss_exp+ss_res, ss_tot]
    [mean(x), var(x), round(mean(y), 2), round(var(y), 3), round(corr(x,y),3), round(p(1:2)',2), round(R2,2)]
end



function p = polynomial_regression(x, y, n);
A = poly_matrix(x, n);
s=warning('off', 'MATLAB:rankDeficientMatrix');
p = A\y;
warning(s);

function mse=evaluate_mse(x, y, p)
yp = polyval(p, x);
mse=sum( (y-yp).^2 );
