function [r, model]=diffusion_1d_residual(model, a, varargin)

f=model.f;
% g=model.g;
K=a(1) * model.K{1};
for i=2:length(a)
    K = K + a(i) * model.K{i};
end
r = f - K*u;
    
% P_I = model.P_I;
% P_B = model.P_B;
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
