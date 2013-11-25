n=300;
[pos,els]=create_mesh_1d(0,1,n);
C=covariance_matrix( pos, @(x1, x2)(matern_covariance(0.1, x1, x2, 0.3, 2)) );
options.correct_var=true;
G_N=mass_matrix( pos, els );
r_i_k=kl_solve_evp( C, G_N, 6, options );
plot( pos, r_i_k)
