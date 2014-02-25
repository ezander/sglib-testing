function model=diffusion_1d_init(varargin)

options=varargin2options(varargin);
[N,options]=get_option(options,'N', 101);
[step_relax,options]=get_option(options,'step_relax', 0.4);
check_unsupported_options(options, mfilename);

num_params = 2;
num_vars = N;

% initialise the model object
model = model_init(num_params, num_vars, ...
    'solve_func', @diffusion_1d_solve, ...
    'step_func', @diffusion_1d_step, ...
    'res_func', @diffusion_1d_residual, ...
    'sol_init_func', @(a)(zeros(N,1)));

model.step_relax = step_relax;

% create mesh and store information in model
[pos,els,bnd]=create_mesh_1d(0, 1, N);
model.pos = pos;
model.els = els;
model.bnd = bnd;

% store the discretised coefficient functions in model
model.a{1}=double(pos<0.5)';
model.a{2}=double(pos>=0.5)';

% store some more FEM stuff in model
%model.u0=zeros(size(pos))';
model.f=ones(size(pos))';
model.g=zeros(size(pos))';
model.K{1}=stiffness_matrix(pos, els, model.a{1});
model.K{2}=stiffness_matrix(pos, els, model.a{2});
[model.P_I,model.P_B]=boundary_projectors( bnd, N );

