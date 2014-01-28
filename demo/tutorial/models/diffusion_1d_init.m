function [state,info]=diffusion_1d_init(varargin)

N=101;
[pos,els,bnd]=create_mesh_1d(0, 1, N);
state.a{1}=double(pos<0.5)';
state.a{2}=double(pos>=0.5)';

state.pos = pos;
state.els = els;
state.bnd = bnd;

state.f=ones(size(pos))';
state.g=zeros(size(pos))';
state.K{1}=stiffness_matrix(pos, els, state.a{1});
state.K{2}=stiffness_matrix(pos, els, state.a{2});

[state.P_I,state.P_B]=boundary_projectors( bnd, N );

info.num_params = 2;
info.num_vars = N;
