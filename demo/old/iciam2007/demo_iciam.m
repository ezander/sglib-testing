%function demo_iciam

disp('Warning: not all parts of this demo work any more, since quite a few ');
disp('         functions were replaced and this demo not updated...' );

addpath( '..' )
%init_demos
clf
clear

% Solving with stochastic operator
global n pos els M %#ok
global p_f m_gam_f m_f lc_f h_f cov_f f_i_alpha I_f%#ok
global p_k m_gam_k m_k lc_k h_k cov_k k_i_alpha I_k%#ok
global p_u m_gam_u I_u %#ok

% loads everything into the global variables
model();


% show the sparsity structure of the deterministic, stochastic-only,
% complete stochastic system matrix
global K_ab K_mu_delta K_mu_iota %#ok
operators();

if iscell(K_ab)
    K_ab_mat=cell2mat(K_ab);
else
    K_ab_mat=K_ab;
end

if 1<0 % don't exectu
    g_i_alpha=0*f_i_alpha;
    [P_B,P_I]=boundary_projectors( [1,n], n );
    [K_ab_mat_s,f_i_alpha_s]=apply_boundary_conditions( K_ab_mat, f_i_alpha(:), g_i_alpha(:), P_B, P_I, size(I_u,1), 'spatial_pos', 2 );
    K_ab_mat=K_ab_mat_s;
    f_i_alpha=reshape(f_i_alpha_s, f_i_alpha );
end

show_sparsity_pattern( K_ab_mat, n );
userwait;

f_beta=compute_pce_rhs( f_i_alpha, I_f, I_u );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% solve
fprintf('\n\n==================================================================================================\n');
%solve_method=6;
solve_method=1;
tol=1e-4; maxiter=50;

m=size(f_beta,2);
[P_B,P_I]=boundary_projectors( [1,n], n );
g_beta=0*f_beta;
[K_ab_mat_x,f_beta_vec]=...
    apply_boundary_conditions_system( K_ab_mat, f_beta(:), g_beta(:), P_B, P_I);

K_ab_mat=K_ab_mat_x;
f_beta=reshape( f_beta_vec, size(f_beta));

u_i_alpha0=K_ab_mat\f_beta(:);
%u_i_alpha0=u_i_alpha0+K_ab_mat\(f_beta(:)-K_ab_mat*u_i_alpha0);
u_i_alpha0=reshape( u_i_alpha0, size(f_beta) );

switch solve_method
    case 1 % use direct solver 
        u_i_alpha=K_ab_mat\f_beta(:);
        itermethod='direct solver'; iter=1; relres=0; flag=0;
    case 2 % pcg with matrices
        M=tkron( K_ab_mat(1:n,1:n), spdiags(multiindex_factorial(I_u),0,size(I_u,1),size(I_u,1)) );
        [u_i_alpha,flag,relres,iter]=pcg( K_ab_mat, f_beta(:), tol, maxiter, M );
        itermethod='pcg(mat)';
    case 3 % use pcg on flat system
        Minv=tkron( spdiags(1./multiindex_factorial(I_u),0,size(I_u,1),size(I_u,1)), inv(K_ab{1,1}) );
        A_func=@(x)(K_ab_mat*x);
        Minv_func=@(x)(Minv*x);
        [u_i_alpha,flag,relres,iter]=pcg( A_func, f_beta(:), tol, maxiter, Minv_func );
        itermethod='pcg(funcs,flat)';
    case 4 % 
        %Minv=tkron( inv(K_mu_delta{1}), inv(sparse(K_mu_delta{2})) );
        Minv_op={ inv(K_mu_delta{1}); inv(sparse(K_mu_delta{2})); {}; {} };
        shape=size(f_beta);
        %Minv_func=@(x)(apply_flat(Minv,x,shape));
        Minv_func=@(x)(apply_flat(Minv_op,x,shape));
        %A_func=@(x)(apply_flat(K_mu_delta,x,shape));
        A_func={@apply_flat,{K_mu_delta,shape},{1,3}};
        [u_i_alpha,flag,relres,iter]=pcg( A_func, f_beta(:), tol, maxiter, Minv_func );
        itermethod='pcg(funcs,mu_delta)';
    case 5
        options.method='gs';
        options.omega=1;
        options.abstol=tol;
        options.reltol=tol;
        options.maxiter=maxiter;
        [u_i_alpha,flag,relres,iter]=solve_stat( K_ab_mat, f_beta(:), options );
        itermethod='mat_decomp(mat)';
    case 6
        options.method='gs';
        options.omega=1;
        options.abstol=tol;
        options.reltol=tol;
        options.maxiter=maxiter;
        options.maxiter=10; % check for divergence(!!)
        
        [u_i_alpha,flag,relres,iter]=solve_mat_decomp_block( K_ab, f_beta, options );
        itermethod='mat_decomp(block,flat)';
    case 7
        options.method='jac';
        options.omega=1;
        options.abstol=tol;
        options.reltol=tol;
        options.trunc_eps=0.0001;

        maxiter=100;
        options.abstol=0.0001;
        options.reltol=0.0001;
        options.trunc_eps=0.000001;
        options.trunc_k=20;
        options.maxiter=maxiter;
        options.relax=1;
        options.algorithm=2;
        
        [U,S,V]=svd(f_beta);
        n=rank(S);
        f1=U(:,1:n)*S(1:n,1:n);
        f2=V(:,1:n);
        F_tp={f1,f2};
        
        %K_mu_delta{1}=1*K_mu_delta{1};
        for i=1:size(K_mu_delta{3},1); K_mu_delta{3}{i}=1.0*K_mu_delta{3}{i}; end
        
        K_mu_delta{2}=sparse(K_mu_delta{2});
        for i=1:size(K_mu_delta{4},1); K_mu_delta{4}{i}=sparse(K_mu_delta{4}{i}); end
        
        [u_i_alpha,flag,relres,iter]=solve_stat_tensor( K_mu_delta, F_tp, options );
        itermethod='mat_decomp(tensor)';
        u_i_alpha=u_i_alpha{1}*u_i_alpha{2}';
    case 99 % statements saved for later use
end


if isvector(u_i_alpha)
    u_i_alpha=reshape( u_i_alpha, size(f_beta) );   
end

fprintf( 'method:   %d\n', solve_method );
fprintf( 'residual: %g\n', normest(operator_apply( K_ab, u_i_alpha(:))-f_beta(:))/normest(f_beta(:)))
fprintf( 'error:    %g\n', normest(u_i_alpha-u_i_alpha0)/normest(u_i_alpha0))
%return

% Show sample realizations and 
show_in_out_samples( pos, k_i_alpha, f_i_alpha, u_i_alpha, I_k, I_f, I_u, 30 );
userwait;
%return

% Show covariances
show_covariances2( pos, k_i_alpha, f_i_alpha, u_i_alpha, I_k, I_f, I_u, cov_k, cov_f );
userwait;

