function [un, model] = electrical_network_picard_iter_step(model, u, p)
% ELECTRICAL_NETWORK_PICARD_ITER_STEP Computes one Picard iteration for the electrical network.
%
% Note: could also be called modified Newton with linearised Jacobian at u=0

% Compute residual at u
[res, model] = electrical_network_residual(model, u, p);

% Compute one Picard iteration step
du = model.A \ res;
un = u + du;
