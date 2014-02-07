function [un, minfo] = electrical_network_picard_iter_step(minfo, u, p)
% ELECTRICAL_NETWORK_PICARD_ITER_STEP Computes one Picard iteration for the electrical network.
%
% Note: could also be called modified Newton with linearised Jacobian at u=0

% Compute residual at u
[res, minfo] = electrical_network_residual(minfo, u, p);

% Compute one Picard iteration step
du = minfo.A \ res;
un = u + du;
