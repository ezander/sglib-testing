function [r, state]=diffusion_1d_residual(state, a, varargin)

f=state.f;
% g=state.g;
K=a(1) * state.K{1};
for i=2:length(a)
    K = K + a(i) * state.K{i};
end
r = f - K*u;
    
% P_I = state.P_I;
% P_B = state.P_B;
% 
% Ki=apply_boundary_conditions_operator( K, P_I );
% fi=apply_boundary_conditions_rhs( K, f, g, P_I, P_B );
% 
% ui=Ki\fi;
% 
% u=apply_boundary_conditions_solution(ui, g, P_I, P_B);
% 
% solve_info=struct();
% solve_info.iter=1;
% solve_info.res=f-K*u;
