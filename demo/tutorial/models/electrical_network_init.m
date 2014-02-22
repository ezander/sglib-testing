function model = electrical_network_init( varargin )
% ELECTRICAL_NETWORK_INIT Initialises the structure that keeps the internal
% model of the electrical network example.
%
% Example:
%    model = electrical_network_init('R', 200, 'f0', 5);

options = varargin2options(varargin);
[R, options] = get_option(options, 'R', 0.01);
[fg, options] = get_option(options, 'fg', 25);
[f0, options] = get_option(options, 'f0', [1; zeros(4,1)]);
[newton, options] = get_option(options, 'newton', true);
check_unsupported_options(options, mfilename);


A=create_network_matrix(R);
num_params = 2;
num_vars = size(A,1);

% initialise the model object
if newton
    step_func = @electrical_network_newton_step;
else
    step_func = @electrical_network_picard_iter_step;
end

model = model_init(num_params, num_vars, ...
    'solve_func', @electrical_network_solve, ...
    'step_func', step_func, ...
    'res_func', @electrical_network_residual, ...
    'sol_init_func', @(a)(zeros(num_vars,1)));

% store electrical network info in model object
model.A = A;
model.f0 = f0;
model.fg = fg;
%model.u0 = zeros(num_vars, 1);



function A=create_network_matrix(R)

% edge - node incidence matrix
B=[0  1  0 -1  0  0;
   1  0  0  0 -1  0; 
   0  0  0 -1  1  0; 
   0  0 -1  1  0  0; 
   1  0 -1  0  0  0;  
   0 -1  1  0  0  0;
  -1  1  0  0  0  0; 
   0  0  1  0 -1  0;
   0  0  0  0  1 -1];

% conductivity of resistors on edges
n = size(B,1);
D = diag((1/R)*ones(1,n));
As = B'*D*B;

% As should be singular
assert(rank(As)<min(size(As)));

% grounding of the last node is achieved by removing the last row and
% column of the matrix
A = As(1:end-1, 1:end-1);
assert(rank(A)==min(size(A)));
