%% Some tests on Richardson's iteration
% See: https://en.wikipedia.org/wiki/Modified_Richardson_iteration

% Setup the matrix to be spd
rand_seed( 29138 );
N=10;
A=rand(N)-0.5;
A=A+A';
b=rand(N,1);
A=A+3*eye(N);

% Compute spectral radius and eigenvalues
clc
lambda=eig(A);
rho=max(abs(lambda));
strvarexpand('Spectral radius of A: $rho$');
lambda_min=min(lambda);
lambda_max=max(lambda);
strvarexpand('Minimal eigenvalue of A: $lambda_min$');
strvarexpand('Maximal eigenvalue of A: $lambda_max$');

% The optimal omega
w_opt=2/(lambda_min+lambda_max);
z_opt=rho*w_opt;

% Look at convergenze for different omega's 
step = 0.02;
N = 50;
for dz=-20*step:step:20*step
    z = z_opt + dz;
    w=z/rho;
    
    x=zeros(size(b));
    for i=1:N
        x=x+w*(b-A*x);
    end
    if abs(dz)<=step/2; fprintf('\nSupposed optimum:\n'); end
    fprintf( 'dz: %g  x:%g  r:%g\n', dz, norm(x)/norm(b\A), norm(b-A*x)/norm(b) );
    if abs(dz)<=step/2; disp(' '); end
end
