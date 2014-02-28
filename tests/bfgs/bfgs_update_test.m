function bfgs_update_test


return

% test 
B2b = B1 + [y, B1*s] * [1/(y'*y), 0; 0, -1/(s'*B1*s)] * [y, B1*s]';
norm(B2-B2b)

U=[y, B1*s];
C=[1/(y'*y), 0; 0, -1/(s'*B1*s)];
B2b = B1 + U*C*U';
norm(B2-B2b)
H2b = H1 - H1 * U * inv(inv(C) + U' * H1 * U) * U' * H1;
norm(H2b*B2-I)
norm(H2*B2-I)

Cinv=[(y'*y), 0; 0, -(s'*B1*s)];
H2b = H1 - H1 * U * inv(Cinv + U' * H1 * U) * U' * H1;
norm(H2b*B2-I)

