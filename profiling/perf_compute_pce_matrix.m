hermite_triple_fast(10);
I_k=multiindex(4,4);
I_f=multiindex(2,3);
[I_k,I_f,I_u]=multiindex_combine( {I_k, I_f}, -1 );
n=10;
k_iota=rand(n,size(I_k,1));
tic;
Delta=compute_pce_matrix( k_iota, I_k, I_u );
toc