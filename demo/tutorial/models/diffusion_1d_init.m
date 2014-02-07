function minfo=diffusion_1d_init(varargin)

options=varargin2options(varargin);
[N,options]=get_option(options,'N', 101);
[step_relax,options]=get_option(options,'step_relax', 0.4);
check_unsupported_options(options, mfilename);

num_params = 2;
num_vars = N;

% initialise the minfo object
minfo = model_init(num_params, num_vars, ...
    'solve_func', @diffusion_1d_solve, ...
    'step_func', @diffusion_1d_step, ...
    'res_func', @diffusion_1d_residual, ...
    'sol_init_func', @(a)(zeros(N,1)));

minfo.step_relax = step_relax;

% create mesh and store information in minfo
[pos,els,bnd]=create_mesh_1d(0, 1, N);
minfo.pos = pos;
minfo.els = els;
minfo.bnd = bnd;

% store the discretised coefficient functions in minfo
minfo.a{1}=double(pos<0.5)';
minfo.a{2}=double(pos>=0.5)';

% store some more FEM stuff in minfo
%minfo.u0=zeros(size(pos))';
minfo.f=ones(size(pos))';
minfo.g=zeros(size(pos))';
minfo.K{1}=stiffness_matrix(pos, els, minfo.a{1});
minfo.K{2}=stiffness_matrix(pos, els, minfo.a{2});
[minfo.P_I,minfo.P_B]=boundary_projectors( bnd, N );

