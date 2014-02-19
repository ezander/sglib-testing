clear
clf

I_u=multiindex(3, 2, 'ordering', 'uqtk');
V_u=gpcbasis_create('p', 'I', I_u);

I_a = multiindex(3, 1, 'ordering', 'uqtk');
V_a=gpcbasis_create('p', 'I', I_a);


L = gpcbasis_triples(V_u, V_u, V_a);


spy2(sum(L(:,:,2:end),3))
