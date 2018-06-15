%% Setup data and data matrix
[x, y] = create_data('line', 50, 'sigma', 1);
plot(x, y, '.')

A = [x ones(size(x))];


%% Solve normal equations (matlab builtin)
q = (A'*A)\(A'*y);
m = q(1);
n = q(2);
plot_line(m, n, 2, 5);


%% Solve directly (matlab builtin)
q = A\y;
plot_line(q(1), q(2), 2, 5);


%% Solve using pseudo-inverse (matlab builtin)
Ap = pinv(A);
q = Ap * y;
plot_line(q(1), q(2), 2, 5);






%% Solve normal equations (Cholesky)
C = A' * A;
L = chol(C)'
spy2(L) % check structure of L
norm(C - L*L') % check decomposition

%% Solve 
% A'*A*q = L*L'*q = A'*y
% => q = L' \ (L \ (A'*y))

z = A' * y;
z = L \ z;  % forward substitution (builtin, automatic)
z = L' \ z; % back substitution (builtin, automatic)

q = z;
plot(x, y, '.')
plot_line(q(1), q(2), 2, 5);



%% Solve normal equations (LDL)
C = A' * A;
[L,D] = ldl(C);
spy2(L) % check structure of L
norm(C - L*D*L') % check decomposition

%% Solve 
% A'*A*q = L*D*L'*q = A'*y
% => q = L' \ (L \ (A'*y))

z = A' * y;
z = L \ z;  % forward substitution (builtin, automatic)
z = D \ z;  % solve with diagonal (z = z ./ diag(D))
z = L' \ z; % back substitution (builtin, automatic)

q = z;
plot(x, y, '.')
plot_line(q(1), q(2), 2, 5);




%% Solve directly using QR 
% Perform QR decomposition and check
[Q,R]=qr(A,0);
spy2(R) % check structure of R
assert_equals(A, Q*R, 'check QR decomposition')
assert_matrix(Q, 'orthogonal', 'Q is orth')  % equivalent: norm(Q'*Q-eye(2))<tol 
assert_matrix(R, 'upper', 'R is right/upper')% equivalent: norm(R-triu(R))==0

% Find q such that min(Q*R*q - y)
% => R'*Q'*Q*R*q = R'*Q'*y
% => R*q = Q'*y
% => q = R \ (Q'*y)

q = R \ (Q'*y);
plot(x, y, '.')
plot_line(q(1), q(2), 2, 5);



%% Solve directly using SVD
% Perform SVD  and check
[U,S,V]=svd(A,0);
assert_equals(A, U*S*V', 'check SVD')
assert_matrix(U, 'orthogonal', 'U is orth')  
assert_matrix(S, 'diagonal', 'S is diagonal')
assert_matrix(V, 'orthogonal', 'V is orth')  

% Find q such that min(U*S*V'*q - y)
% => V*S*U' *U*S*V'*q = V*S*U'*y
% => S*V'q = U'*y
% => q = S \ (U'*y)

q = V * (S \ (U'*y));
plot(x, y, '.')
plot_line(q(1), q(2), 2, 5);



%% Solve directly using SVD and pseudo inverse
% Perform SVD  and check
[U,S,V]=svd(A,0);
assert_equals(A, U*S*V', 'check SVD')
assert_matrix(U, 'orthogonal', 'U is orth')  
assert_matrix(S, 'diagonal', 'S is diagonal')
assert_matrix(V, 'orthogonal', 'V is orth')  

% Find q such that min(U*S*V'*q - y)
% => V*S*U' *U*S*V'*q = V*S*U'*y
% => S*V'q = U'*y
% => q = S \ (U'*y)

q = V * (S \ (U'*y));
plot(x, y, '.')
plot_line(q(1), q(2), 2, 5);







