N = 10000;
x = linspace(0, 1, N);
A = poly_matrix(x, 5, 1)/sqrt(N);
(A' * A) - hilb(6)

cond(A)
cond(hilb(6))^0.5

