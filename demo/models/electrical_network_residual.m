function [res, model] = electrical_network_residual(model, u, p)
% ELECTRICAL_NETWORK_RESIDUAL function to compute the residuum
% (preconditioned).
%
% Each '*_residual' function has the system model as first parameter, the
% approximate solution, for which the residual is computed, as the second
% parameter, and the vector of variable or stochstic parameters as the
% third parameter.

A = model.A;
fg = model.fg;
f0 = model.f0;

% Compute the residual of the electrical network
res = (fg + p(2)).*f0 - (A*u + (p(1)+2)*(u'*u)*u);
