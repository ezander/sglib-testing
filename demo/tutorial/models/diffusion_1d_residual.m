function [r, minfo]=diffusion_1d_residual(minfo, a, varargin)

f=minfo.f;
% g=minfo.g;
K=a(1) * minfo.K{1};
for i=2:length(a)
    K = K + a(i) * minfo.K{i};
end
r = f - K*u;
    
% P_I = minfo.P_I;
% P_B = minfo.P_B;
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
