function Delta=legendre_triples(n, l)
% LEGENDRE_TRIPLES Return triples E(P_iP_jP_k) for i,j<n and k<l
V_Y=gpcbasis_create('P', 'm', 1, 'p', n-1);
V_X=gpcbasis_create('P', 'm', 1, 'p', l-1);
Delta = gpcbasis_triples(V_Y, V_Y, V_X);


