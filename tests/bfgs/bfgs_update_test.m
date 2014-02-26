function bfgs_update_test

%
N=6;
A1 = rand(N);
B1 = A1'*A1+eye(N);
%B1 = A1;
B1 = A1'*A1-0.5*eye(N);

H1 = inv(B1);
I = eye(N);

y = rand(N,1);
s = rand(N,1);

clc
unittest_updates(B1, H1, y, s);

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

function unittest_updates(B1, H1, y, s)

[B2, H2] = update_dfp(B1, H1, y, s);
assert_matrix(B2*H2, 'identity', 'id_dfp');

[B2, H2] = update_bsgs(B1, H1, y, s);
assert_matrix(B2*H2, 'identity', 'id_bfgs');

[B2, H2] = update_sr1(B1, H1, y, s);
assert_matrix(B2*H2, 'identity', 'id_sr1');

[B2, H2] = update_broyden(B1, H1, y, s);
assert_matrix(B2*H2, 'identity', 'id_broyden');


function [B2, H2] = update_dfp(B1, H1, y, s)
I=eye(size(B1));
rho = 1/(y'*s);
B2 = (I-rho*y*s')*B1*(I-rho*s*y') + rho*(y*y');
H2 = H1 + rho*(s*s') - (H1*y)*(H1*y)'/(y'*H1*y);

function [B2, H2] = update_bsgs(B1, H1, y, s)
I=eye(size(B1));
rho = 1/(y'*s);
B2 = B1 + rho*(y*y') - (B1*s)*(B1*s)'/(s'*B1*s);
H2 = (I-rho*s*y')*H1*(I - rho*y*s') + rho*(s*s');

function [B2, H2] = update_sr1(B1, H1, y, s)
B2 = B1 + (y-B1*s)*(y-B1*s)'/((y-B1*s)'*s);
H2 = H1 + (s-H1*y)*(s-H1*y)'/((s-H1*y)'*y);

function [B2, H2] = update_broyden(B1, H1, y, s)
B2 = B1 + (y-B1*s)/(s'*s)*s';
H2 = H1 + (s-H1*y)*s'*H1/(s'*H1*y);


