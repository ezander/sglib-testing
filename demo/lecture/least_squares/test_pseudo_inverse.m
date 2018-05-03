function test_pseudo_inverse


%% Setup data and matrix
[x, ~] = create_data('line', 105, 'sigma', 1);
A = poly_matrix(x, 5);
% The following line makes A rank deficient and thus makes it more difficult for the algorithm
%A = [A, A] % comment out to have full row rank

% check properties (pinv)
Ap = pinv(A);
check_pinv_properties(A, Ap)

%% compute via QR
[Q,R] = qr(A, 0);
Ap2 = R\Q';
norm(Ap2 - Ap)
check_pinv_properties(A, Ap2)

%% compute via QR
[Q,R] = qr(A, 0);
Ap2 = (R'*R)\A'
norm(Ap2 - Ap)
check_pinv_properties(A, Ap2)

%% compute via SVD
[U,S,V] = svd(A, 0);
Sp = 1./S;
Sp(Sp>1e10)=0; % very crude (better/more general: > 1/ (eps * max(size(A)) * max(S(:))) )
Ap3 = V*Sp*U';
norm(Ap3 - Ap)
check_pinv_properties(A, Ap3)



