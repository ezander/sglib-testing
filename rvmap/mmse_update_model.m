function [xn_i_beta, V_xn]=mmse_update_model(x_func, y_func, V_xy, ym_beta, V_ym, p_phi, p_int_mmse, p_xn, p_int_proj)
% MMSE_UPDATE_MODEL Update a stochastic model by the MMSE method.
%   [XN_I_BETA, V_XN]=MMSE_UPDATE_MODEL(X_FUNC, Y_FUNC, V_XY, YM_BETA,
%   V_YM, P_PHI, P_INT_MMSE, P_XN, P_INT_PROJ) computes the "updated" model
%   for X. X_FUNC represents the random variable X on input, defined on the
%   germ given by V_XY. Y represents the measurement random variable, i.e.
%   the measurement operator applied to the state of the system (usually
%   the composition of some M_FUNC representing the measurment op and some
%   U_FUNC representing the system state random variable, which must be
%   defined on the same germ V_XY as X_FUNC). YM_BETA and V_YM define a GPC
%   of the actual measurements. Usually, the mean of YM_BETA will be what
%   was actually measured and the GPC will model this measurement plus some
%   error model. 
%   P_PHI is the degree of the polynomials used to compute the MMSE
%   estimator, P_INT_MMSE is the order of integration used there in. 
%   MMSE_UPDATE_MODEL Long description of mmse_update_model. P_XN is the
%   degree of the polynomials in the GPC of the updated model XN.
%   The GPC of the updated model is returned in XN_I_BETA and V_XN.
%
% Example (<a href="matlab:run_example mmse_update_model">run</a>)
%
% See also

%   Elmar Zander
%   Copyright 2014, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

warning('This function is obsolete. Don''t use any more!!!');

[phi_j_delta,V_phi]=mmse_estimate(x_func, y_func, V_xy, p_phi, p_int_mmse);

phi_func = gpc_function(phi_j_delta, V_phi);
ym_func = gpc_function(ym_beta, V_ym);
xn_func = funcompose(ym_func, phi_func);

V_xn = gpcbasis_create(V_ym, 'p', p_xn);
xn_i_beta = gpc_projection(xn_func, V_xn, p_int_proj);
