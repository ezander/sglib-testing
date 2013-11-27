n=7
[xc, wc]=clenshaw_curtis_rule(n+1);

format short g
format compact

[xn, wn]=newton_cotes_rule(n);
xc2=cos(0.5*(1-xn)*pi);

wc2=interpolatory_weights(xc2);

xc
xc2
wc'
wc2'



V_x=gpcbasis_create('L', 'm', 1, 'p', 3);
x_i_alpha = [2 -4 2 0]
gpc_evaluate(x_i_alpha, V_x, [1,2,3,4])

V_y=gpcbasis_create('L', 'm', 1, 'p', 3);
y_j_beta = [6 -18 18 -6];
gpc_evaluate(y_j_beta, V_y, [1,2,3,4])

%%
A=polysys_rc2coeffs(polysys_recur_coeff('L', 3))'
c=A\[1;0;0;0]
A*c
