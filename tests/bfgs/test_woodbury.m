function test_woodbury


N=6;
k=3;
A = rand(N);
U = rand(N,k);
V = rand(N,k);
X = rand(N,2);

A2 = A + U*V';
Ainv = operator_from_matrix_solve(A);

clc
format compact; format short g;
A2\X
woodbury_solve(Ainv, U, V, X)





function Y=woodbury_solve(Ainv, U, V, X)
k = size(U,2);
M = eye(k) + V' * operator_apply(Ainv, U);
Z1 = operator_apply(Ainv, X);
Z2 = operator_apply(Ainv, U * (M \ (V' * Z1)));
Y = Z1 - Z2;
