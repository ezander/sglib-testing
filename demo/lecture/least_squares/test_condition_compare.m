function test_condition_compare
% compare using normal equations and QR approach


%% Setup data and data matrix
[x, y] = create_data('line', 105, 'sigma', 1);
plot(x, y, '.')
A = poly_matrix(x, 11) % change poly degree here and see how solution changes


%% Solve using QR
cond(A)
[Q,R]=qr(A,0);
q = R\(Q'*y);
plot_poly_curve(q, min(x), max(x), 'color', 'g')



%% Solve Normal equations
cond(A'*A)
[L,D]=ldl(A'*A);
q = L'\(D\(L\(A'*y)));
%q(end) = q(end)+0.1;
plot_poly_curve(q, min(x), max(x), 'color', 'r')

