clear
I_u=multiindex(4, 3, 'ordering', 'uqtk');
V_u=gpcbasis_create('p', 'I', I_u);

I_a = multiindex(4, 1, 'ordering', 'uqtk');
V_a=gpcbasis_create('p', 'I', I_a);


L = gpc_triples(V_u, V_u, V_a);


subplot(2,3,1); spy2(sum(L(:,:,2:end),3))
for i=1:5
    subplot(2,3,i+1); spy2(L(:,:,i))
end
