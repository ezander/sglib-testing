function [un, minfo] = electrical_network_newton_step(minfo, u, p)
% ELECTRICAL_NETWORK_NEWTON_STEP Computes one Newton step for the electrical network.

% Compute residual at u
[res, minfo] = electrical_network_residual(minfo, u, p);

% Compute Jacobian at u
A = minfo.A;
I = eye(size(A));
J = -(A + (p(1)+2) * (2*(u*u') + (u'*u)*I));

%du = 0.000001 * rand(length(u),1);
%res2 = electrical_network_residual(minfo, u+du, p);
%norm((res2 - res) - J * du)/norm(du)

% Compute Newton step
du = -J \ res;
un = u + du;