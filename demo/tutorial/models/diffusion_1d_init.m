function state=diffusion_1d_init(varargin)

options=varargin2options(varargin);
[N,options]=get_option(options,'N', 101);
[step_relax,options]=get_option(options,'step_relax', 0.4);
check_unsupported_options(options, mfilename);

num_params = 2;
num_vars = N;

% initialise the state object
state = model_init(num_params, num_vars, @diffusion_1d_solve, @diffusion_1d_step);

state.step_relax = step_relax;

% create mesh and store information in state
[pos,els,bnd]=create_mesh_1d(0, 1, N);
state.pos = pos;
state.els = els;
state.bnd = bnd;

% store the discretised coefficient functions in state
state.a{1}=double(pos<0.5)';
state.a{2}=double(pos>=0.5)';

% store some more FEM stuff in state
state.u0=zeros(size(pos))';
state.f=ones(size(pos))';
state.g=zeros(size(pos))';
state.K{1}=stiffness_matrix(pos, els, state.a{1});
state.K{2}=stiffness_matrix(pos, els, state.a{2});
[state.P_I,state.P_B]=boundary_projectors( bnd, N );

