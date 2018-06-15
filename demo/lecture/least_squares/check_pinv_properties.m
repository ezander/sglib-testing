function check_pinv_properties(A, Ap)
[norm(A*(Ap*A)-A),
norm((Ap*A)*Ap-Ap),
%norm((Ap*A)*A'-A')
norm((Ap*A)-(Ap*A)'),
norm((A*Ap)-(A*Ap)')]
