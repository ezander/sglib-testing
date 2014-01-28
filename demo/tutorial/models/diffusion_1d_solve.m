function [u, a, pos]=diffusion_1d_solve(a1, a2)

N=101;
[pos,els,bnd]=create_mesh_1d(0, 1, N);

a=(a1*double(pos<0.5) + a2*double(pos>=0.5))';
f=ones(size(a));
g=zeros(size(a));
K=stiffness_matrix(pos, els, a);

[P_I,P_B]=boundary_projectors( bnd, N );
Ki=apply_boundary_conditions_operator( K, P_I );
fi=apply_boundary_conditions_rhs( K, f, g, P_I, P_B );

ui=Ki\fi;

u=apply_boundary_conditions_solution(ui, g, P_I, P_B);
