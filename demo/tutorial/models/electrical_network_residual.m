function [res, minfo] = electrical_network_residual(minfo, u, p)
% ELECTRICAL_NETWORK_RESIDUAL function to compute the residuum
% (preconditioned).
%
% Each '*_residual' function has the system minfo as first parameter, the
% approximate solution, for which the residual is computed, as the second
% parameter, and the vector of variable or stochstic parameters as the
% third parameter.

A = minfo.A;
fg = minfo.fg;
f0 = minfo.f0;

% Compute the residual of the electrical network
res = (fg + p(2)).*f0 - (A*u + (p(1)+2)*(u'*u)*u);
